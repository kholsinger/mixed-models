data {
  int<lower=0> n_species;
  int<lower=0> n_sites;
  int<lower=0> n_obs;
  int<lower=0> species[n_obs];
  int<lower=0> site[n_obs];
  real<lower=0> intercept_scale;
  real<lower=0> beta_scale;
  real<lower=0> sigma_scale;
  vector[n_obs] x;
  vector[n_obs] y;
}

parameters {
  real beta_0;
  real beta_1;
  real<lower=0> sigma;
  real<lower=0> sigma_species;
  real<lower=0> sigma_species_site;
  vector[n_species] beta_species;
  matrix[n_species,n_sites] beta_species_site;
}

transformed parameters {
  vector[n_obs] mu_obs;
  vector[n_species] mu_species;
  matrix[n_species,n_sites] mu_species_site;

  for (j in 1:n_species) {
    mu_species[j] = beta_0 + beta_species[j];
    for (k in 1:n_sites) {
      mu_species_site[j,k] = mu_species[j] + beta_species_site[j,k];
    }
  }
  for (i in 1:n_obs) {
    mu_obs[i] = mu_species_site[species[i], site[i]]
                + beta_1*x[i];
  }
}

model {
  y ~ normal(mu_obs, sigma);
  beta_0 ~ normal(0.0, intercept_scale);
  beta_species ~ normal(0.0, sigma_species);
  to_vector(beta_species_site) ~ normal(0.0, sigma_species_site);
  beta_1 ~ normal(0.0, beta_scale);
  sigma ~ cauchy(0.0, sigma_scale);
  sigma_species ~ gamma(1.0, 1.0); //cauchy(0.0, sigma_scale);
  sigma_species_site ~ gamma(1.0, 1.0); //cauchy(0.0, sigma_scale);
}
