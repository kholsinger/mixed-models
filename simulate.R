library(rstanarm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

rm(list=ls())

n_sims <- 100

n_species <- 5
n_sites <- 3
n_indiv <- 20

beta_1 <- 0.5
sigma <- 0.25

sigma_species <- 1.0
sigma_species_site <- 0.75

sigma_x <- 1.0

intercept_scale <- 1.0
beta_scale <- 1.0
sigma_scale <- 5.0
regularize <- 1.0

stan_pars <- c("beta_0",
               "beta_1",
               "sigma",
               "sigma_species",
               "sigma_species_site")

cover <- function(x, target) {
  interval <- quantile(x, c(0.1, 0.9))
  if ((interval[1] > target) || (interval[2] < target)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

bias <- function(parm, results) {
  dat <- subset(results, parameter==parm)
  res <- paste(round(mean(dat$mean_lmer - dat$value), 3), "(lmer) ",
               round(mean(dat$mean_stan - dat$value), 3), "(stan)")
  return(res)
}

rmse <- function(parm, results) {
  dat <- subset(results, parameter==parm)
  res <- paste(round(sqrt(mean(dat$mean_lmer - dat$value)^2), 3), "(lmer) ",
               round(sqrt(mean(dat$mean_stan - dat$value)^2), 3), "(stan)")
  return(res)
}

coverage <- function(parm, results) {
  dat <- subset(results, parameter==parm)
  res <- paste(round(sum(dat$cover_lmer)/nrow(dat), 3), "(lmer) ",
               round(sum(dat$cover_stan)/nrow(dat), 3), "(stan)")
  return(res)
}

stan_pars <- c("beta_0",
               "beta_1",
               "sigma",
               "sigma_species",
               "sigma_species_site")

cover <- function(x, target) {
  interval <- quantile(x, c(0.1, 0.9))
  if ((interval[1] > target) || (interval[2] < target)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

bias <- function(parameter, results) {
  dat <- subset(results, parameter==parameter)
  res <- paste(round(mean(dat$mean_lmer - dat$value), 3), "(lmer) ",
               round(mean(dat$mean_stan - dat$value), 3), "(stan)")
  return(res)
}

rmse <- function(parameter, results) {
  dat <- subset(results, parameter==parameter)
  res <- paste(round(sqrt(mean(dat$mean_lmer - dat$value)^2), 3), "(lmer) ",
               round(sqrt(mean(dat$mean_stan - dat$value)^2), 3), "(stan)")
  return(res)
}

coverage <- function(parameter, results) {
  dat <- subset(results, parameter==parameter)
  res <- paste(round(sum(dat$cover_lmer)/nrow(dat), 3), "(lmer) ",
               round(sum(dat$cover_stan)/nrow(dat), 3), "(stan)")
  return(res)
}

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

iterate_comparison <- function(n_species,
                               n_sites,
                               n_indiv,
                               beta_1,
                               sigma,
                               sigma_species,
                               sigma_species_site)
{
  dat <- generate_data(n_species,
                       n_sites,
                       n_indiv,
                       beta_1,
                       sigma,
                       sigma_species,
                       sigma_species_site)

  fit.lmer <- stan_lmer(y ~ x + (1|species/site),
                        data=dat,
                        prior=normal(0.0, beta_scale),
                        prior_intercept=normal(0.0, intercept_scale),
                        prior_aux=cauchy(0.0, sigma_scale),
                        prior_covariance=decov(regularization=regularize,
                                               concentration=1.0,
                                               shape=1.0,
                                               scale=1.0),
                        adapt_delta=0.99)

  stan_data <- list(species=as.numeric(dat$species),
                    site=as.numeric(dat$site),
                    x=dat$x,
                    y=dat$y,
                    n_species=n_species,
                    n_sites=n_sites,
                    n_obs=n_species*n_sites*n_indiv,
                    intercept_scale=intercept_scale,
                    beta_scale=beta_scale,
                    sigma_scale=sigma_scale)
  fit.stan <- stan(file="simulate.stan",
                   data=stan_data,
                   pars=stan_pars,
                   control=list(adapt_delta=0.99,
                                max_treedepth=20))

  stan_lmer_pars <- c("(Intercept)",
                      "x",
                      "sigma",
                      "Sigma[site:species:(Intercept),(Intercept)]",
                      "Sigma[species:(Intercept),(Intercept)]")

  fit.lmer.stats <- as.data.frame(fit.lmer)[, stan_lmer_pars]
  fit.stan.stats <- extract(fit.stan, pars=stan_pars)

  beta_0_mean <- c(mean(fit.lmer.stats$`(Intercept)`),
                   mean(fit.stan.stats$beta_0))
  beta_0_cover <- c(cover(fit.lmer.stats$`(Intercept)`, 0.0),
                    cover(fit.stan.stats$beta_0, 0.0))
  beta_1_mean <- c(mean(fit.lmer.stats$x),
                   mean(fit.stan.stats$beta_1))
  beta_1_cover <- c(cover(fit.lmer.stats$x, beta_1),
                    cover(fit.stan.stats$beta_1, beta_1))
  sigma_mean <- c(mean(fit.lmer.stats$sigma),
                  mean(fit.stan.stats$sigma))
  sigma_cover <- c(cover(fit.lmer.stats$sigma, sigma),
                   cover(fit.stan.stats$sigma, sigma))
  sigma_species_mean <- c(mean(fit.lmer.stats$`Sigma[species:(Intercept),(Intercept)]`),
                          mean(fit.stan.stats$sigma_species))
  sigma_species_cover <- c(cover(fit.lmer.stats$`Sigma[species:(Intercept),(Intercept)]`, sigma_species),
                           cover(fit.stan.stats$sigma_species, sigma_species))
  sigma_species_site_mean <- c(mean(fit.lmer.stats$`Sigma[site:species:(Intercept),(Intercept)]`),
                               mean(fit.stan.stats$sigma_species_site))
  sigma_species_site_cover <- c(cover(fit.lmer.stats$`Sigma[site:species:(Intercept),(Intercept)]`, sigma_species_site),
                                cover(fit.stan.stats$sigma_species_site, sigma_species_site))
  return(list(beta_0_mean=beta_0_mean,
              beta_0_cover=beta_0_cover,
              beta_1_mean=beta_1_mean,
              beta_1_cover=beta_1_cover,
              sigma_mean=sigma_mean,
              sigma_cover=sigma_cover,
              sigma_species_mean=sigma_species_mean,
              sigma_species_cover=sigma_species_cover,
              sigma_species_site_mean=sigma_species_site_mean,
              sigma_species_site_cover=sigma_species_site_cover))
}

n_pars <- length(stan_pars)
parameter <- character(n_sims*n_pars)
value <- numeric(n_sims*n_pars)
mean_lmer <- numeric(n_sims)
mean_stan <- numeric(n_sims)
cover_lmer <- numeric(n_sims)
cover_stan <- numeric(n_sims)
for (i in 1:n_sims) {
  sink("status.txt")
  cat("Iteration: ", i)
  sink()
  res <- iterate_comparison(n_species,
                            n_sites,
                            n_indiv,
                            beta_1,
                            sigma,
                            sigma_species,
                            sigma_species_site)
  parameter[(i-1)*n_pars + 1] <- "beta_0"
  value[(i-1)*n_pars + 1] <- 0.0
  mean_lmer[(i-1)*n_pars + 1] <- res$beta_0_mean[1]
  mean_stan[(i-1)*n_pars + 1] <- res$beta_0_mean[2]
  cover_lmer[(i-1)*n_pars + 1] <- res$beta_0_cover[1]
  cover_stan[(i-1)*n_pars + 1] <- res$beta_0_cover[2]
  parameter[(i-1)*n_pars + 2] <- "beta_1"
  value[(i-1)*n_pars + 2] <- beta_1
  mean_lmer[(i-1)*n_pars + 2] <- res$beta_1_mean[1]
  mean_stan[(i-1)*n_pars + 2] <- res$beta_1_mean[2]
  cover_lmer[(i-1)*n_pars + 2] <- res$beta_1_cover[1]
  cover_stan[(i-1)*n_pars + 2] <- res$beta_1_cover[2]
  parameter[(i-1)*n_pars + 3] <- "sigma"
  value[(i-1)*n_pars + 3] <- sigma
  mean_lmer[(i-1)*n_pars + 3] <- res$sigma_mean[1]
  mean_stan[(i-1)*n_pars + 3] <- res$sigma_mean[2]
  cover_lmer[(i-1)*n_pars + 3] <- res$sigma_cover[1]
  cover_stan[(i-1)*n_pars + 3] <- res$sigma_cover[2]
  parameter[(i-1)*n_pars + 4] <- "sigma_species"
  value[(i-1)*n_pars + 4] <- sigma_species
  mean_lmer[(i-1)*n_pars + 4] <- res$sigma_species_mean[1]
  mean_stan[(i-1)*n_pars + 4] <- res$sigma_species_mean[2]
  cover_lmer[(i-1)*n_pars + 4] <- res$sigma_species_cover[1]
  cover_stan[(i-1)*n_pars + 4] <- res$sigma_species_cover[2]
  parameter[(i-1)*n_pars + 5] <- "sigma_species_site"
  value[(i-1)*n_pars + 5] <- sigma_species_site
  mean_lmer[(i-1)*n_pars + 5] <- res$sigma_species_site_mean[1]
  mean_stan[(i-1)*n_pars + 5] <- res$sigma_species_site_mean[2]
  cover_lmer[(i-1)*n_pars + 5] <- res$sigma_species_site_cover[1]
  cover_stan[(i-1)*n_pars + 5] <- res$sigma_species_site_cover[2]
}

results <- data.frame(parameter=parameter,
                      value=value,
                      mean_lmer=mean_lmer,
                      mean_stan=mean_stan,
                      cover_lmer=cover_lmer,
                      cover_stan=cover_stan)

sink("results.txt", append=TRUE)
cat("n_sims:    ", n_sims, "\n",
    "n_species: ", n_species, "\n",
    "n_sites:   ", n_sites, "\n",
    "n_indiv:   ", n_indiv, "\n\n", sep="")
cat("intercept_scale ", intercept_scale, "\n",
    "beta_scale ", beta_scale, "\n",
    "sigma_scale ", sigma_scale, "\n",
    "regularize ", regularize, "\n\n", sep="")
cat("cauchy(0.0, sigma_scale) for random effects in stan",
    "vectorized version\n\n")
for (parm in stan_pars) {
  cat(parm, "\n",
      "      bias: ", bias(parm, results), "\n",
      "      rmse: ", rmse(parm, results), "\n",
      "  coverage: ", coverage(parm, results), "\n")
}
sink()
