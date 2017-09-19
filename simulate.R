library(rstanarm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

rm(list=ls())

n_species <- 5
n_sites <- 3
n_indiv <- 20

beta_1 <- 0.5
sigma <- 0.25

sigma_species <- 1.0
sigma_species_site <- 0.75

sigma_x <- 1.0

beta_scale <- 1.0
sigma_scale <- 5.0

generate_data <- function(n_species,
                          n_sites,
                          n_indiv,
                          beta_1,
                          sigma,
                          sigma_species,
                          sigma_species_site)
{
  species <- numeric(0)
  site <- numeric(0)
  x <- numeric(0)
  y <- numeric(0)
  beta_0_species <- rnorm(n_species, 0.0, sigma_species)
  for (i in 1:n_species) {
    species <- c(species, rep(i, n_sites*n_indiv))
    beta_0_species_site <- rnorm(n_sites,
                                 beta_0_species[i],
                                 sigma_species_site)
    for (j in 1:n_sites) {
      site <- c(site, rep(j, n_indiv))
      x_tmp <-  rnorm(n_indiv, 0.0, sigma_x)
      x <- c(x, x_tmp)
      mu <- beta_0_species_site[j] + beta_1*x_tmp
      for (k in 1:n_indiv) {
        y <- c(y, rnorm(1, mu[k], sigma))
      }
    }
  }
  dat <- data.frame(species=as.factor(species),
                    site=as.factor(site),
                    x=x,
                    y=y)
  return(dat)
}

dat <- generate_data(n_species,
                     n_sites,
                     n_indiv,
                     beta_1,
                     sigma,
                     sigma_species,
                     sigma_species_site)

fit.lmer <- stan_lmer(y ~ x + (1|species/site), data=dat, adapt_delta=0.99)

stan_data <- list(species=as.numeric(dat$species),
                  site=as.numeric(dat$site),
                  x=dat$x,
                  y=dat$y,
                  n_species=n_species,
                  n_sites=n_sites,
                  n_obs=n_species*n_sites*n_indiv,
                  beta_scale=beta_scale,
                  sigma_scale=sigma_scale)
stan_pars <- c("beta_0",
               "beta_1",
               "sigma",
               "sigma_species",
               "sigma_species_site")
fit.stan <- stan(file="simulate.stan",
                 data=stan_data,
                 pars=stan_pars,
                 control=list(adapt_delta=0.99))

print(summary(fit.lmer,
              pars=c("(Intercept)",
                     "x",
                     "sigma",
                     "Sigma[site:species:(Intercept),(Intercept)]",
                     "Sigma[species:(Intercept),(Intercept)]"),
              digits=3))
print(fit.stan, digits=3)

