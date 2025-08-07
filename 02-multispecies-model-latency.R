#Modelling species-specific differences in latency times

set.seed(1234)
source("00-data-prep-exploration.R")

{
  #Model fitting
  # Hurdle gamma with zeroes in function of the species, full model
  
  # Define priors
  {
    priors <- c(
      # Prior for the intercept (main model)
      prior(normal(0, 2), class = "Intercept"),
      # Priors for main fixed effects (default class = "b")
      prior(normal(0, 2), class = "b"),
      # Prior for shape parameter (Gamma)
      prior(gamma(0.01, 0.01), class = "shape"),
      # Prior for group-level SD (random effect of area)
      prior(exponential(1), class = "sd"),
      
      # Prior for the intercept (hurdle part)
      prior(normal(0, 1), class = "Intercept", dpar = "hu"),
      # Priors for fixed effects in the hurdle part
      prior(normal(0, 1), class = "b", dpar = "hu"),
      # Random effect SD for hu model
      prior(exponential(1), class = "sd", dpar = "hu")
    )
  }
  
  # Fit full model
  {
    mod.a.latency <- brm(bf(latency.time.adjusted ~ elevation.std*sp + shannon.std*sp + sex +
                              (1 | area),
                            hu ~ elevation.std + sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    # save(mod.a.latency, file = "output/mod.a.latency.RData", compress = F)
  }
  
  # Multicollinearity check
  {
    library(posterior)
    # Extract posterior draws as a data frame
    post <- as_draws_df(mod.a.latency)
    names(post)
    cor(post$b_elevation.std, post$b_shannon.std) # if <= 0.3 it's ok
    plot(post$b_elevation.std, post$b_shannon.std,
         xlab = "Slope: elevation.std",
         ylab = "Slope: shannon.std",
         main = "Posterior correlation between elevation and Shannon effects")
    posterior_cor <- cor(post[, grep("^b_", names(post))])
    round(posterior_cor, 2)
    # View elevation-related terms vs shannon-related terms
    elevation_terms <- grep("elevation.std", colnames(posterior_cor), value = TRUE)
    shannon_terms   <- grep("shannon.std", colnames(posterior_cor), value = TRUE)
    # Print correlations between elevation and shannon slopes
    posterior_cor[elevation_terms, shannon_terms]
  }
  
  # Sensitivity analyses for slope prior
  {
    priors.sens <- c(
      # Prior for the intercept (main model)
      prior(normal(0, 2), class = "Intercept"),
      # Priors for main fixed effects (default class = "b")
      prior(normal(0, 5), class = "b"),
      # Prior for shape parameter (Gamma)
      prior(gamma(0.01, 0.01), class = "shape"),
      # Prior for group-level SD (random effect of area)
      prior(exponential(1), class = "sd"),
      
      # Prior for the intercept (hurdle part)
      prior(normal(0, 1), class = "Intercept", dpar = "hu"),
      # Priors for fixed effects in the hurdle part
      prior(normal(0, 5), class = "b", dpar = "hu"),
      # Random effect SD for hu model
      prior(exponential(1), class = "sd", dpar = "hu")
      
    )
    mod.a.latency.sens <- brm(bf(latency.time.adjusted ~ elevation.std*sp + shannon.std*sp + sex +
                                   (1 | area),
                                 hu ~ elevation.std + sp + (1 | area)),
                              family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                              data=d,
                              chains=4, cores=4,
                              iter=5000, warmup=2000,
                              prior = priors.sens,
                              # thin = 10,
                              control = list(adapt_delta = 0.9999, max_treedepth = 15),
                              seed = 1234)
    # save(mod.a.latency.sens, file = "output/mod.a.latency.sens.RData", compress = F)
  }
  
  # print results
  {
    print(summary(mod.a.latency))
    print(summary(mod.a.latency.sens))
    print(posterior_cor[elevation_terms, shannon_terms])
  }
}

## ---- Model selection ----
{
  # Comparing models through leave-one-out cross validation
  
  # Fitting nested models
  {
    mod.b.latency <- brm(bf(latency.time.adjusted ~ elevation.std*sp +
                              (1 | area),
                            hu ~ elevation + sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.c.latency <- brm(bf(latency.time.adjusted ~ elevation.std*sp +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.d.latency <- brm(bf(latency.time.adjusted ~ elevation.std + sp +
                              (1 | area),
                            hu ~ elevation + sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.e.latency <- brm(bf(latency.time.adjusted ~ elevation.std + sp +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.f.latency <- brm(bf(latency.time.adjusted ~ sp +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.g.latency <- brm(bf(latency.time.adjusted ~ 1 +
                              (1 | area),
                            hu ~ 1 + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         # prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.h.latency <- brm(bf(latency.time.adjusted ~ elevation.std * shannon.std + sp + sex +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.i.latency <- brm(bf(latency.time.adjusted ~ elevation.std * shannon.std * sp + sex +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    mod.j.latency <- brm(bf(latency.time.adjusted ~ elevation.std * shannon.std +
                              (1 | area),
                            hu ~ sp + (1 | area)),
                         family=hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
                         data=d,
                         chains=4, cores=4,
                         iter=5000, warmup=2000,
                         prior = priors,
                         # thin = 10,
                         control = list(adapt_delta = 0.9999, max_treedepth = 15),
                         seed = 1234)
    
    
    # save(mod.b.latency, file = "output/mod.b.latency.RData", compress = F)
    # save(mod.c.latency, file = "output/mod.c.latency.RData", compress = F)
    # save(mod.d.latency, file = "output/mod.d.latency.RData", compress = F)
    # save(mod.e.latency, file = "output/mod.e.latency.RData", compress = F)
    # save(mod.f.latency, file = "output/mod.f.latency.RData", compress = F)
    # save(mod.g.latency, file = "output/mod.g.latency.RData", compress = F)
    # save(mod.h.latency, file = "output/mod.h.latency.RData", compress = F)
    # save(mod.i.latency, file = "output/mod.i.latency.RData", compress = F)
    # save(mod.j.latency, file = "output/mod.j.latency.RData", compress = F)
    
  }
  
  # LOO cross-validation
  # load("output/mod.a.latency.RData")
  # load("output/mod.b.latency.RData")
  # load("output/mod.c.latency.RData")
  # load("output/mod.d.latency.RData")
  # load("output/mod.e.latency.RData")
  # load("output/mod.f.latency.RData")
  # load("output/mod.g.latency.RData")
  # load("output/mod.h.latency.RData")
  # load("output/mod.i.latency.RData")
  # load("output/mod.j.latency.RData")
  {
    loo_a_full <- loo(mod.a.latency, reloo=TRUE, seed = 1234) # Use it only if k > 0.7 warnings in loo() output
    loo_b <- loo(mod.b.latency, seed = 1234)
    loo_c <- loo(mod.c.latency, seed = 1234)
    loo_d <- loo(mod.d.latency, seed = 1234)
    loo_e <- loo(mod.e.latency, seed = 1234)
    loo_f <- loo(mod.f.latency, seed = 1234)
    loo_g <- loo(mod.g.latency, seed = 1234)
    loo_h <- loo(mod.h.latency, seed = 1234)
    loo_i <- loo(mod.i.latency, reloo = TRUE, seed = 1234)
    loo_j <- loo(mod.j.latency, reloo = TRUE, seed = 1234)
  }
  
  # print results
  {
    print("ELPD mod.a")
    print(loo_a_full)
    print("ELPD mod.b")
    print(loo_b)
    print("ELPD mod.c")
    print(loo_c)
    print("ELPD mod.d")
    print(loo_d)
    print("ELPD mod.e")
    print(loo_e)
    print("ELPD mod.f")
    print(loo_f)
    print("ELPD mod.g")
    print(loo_g)
    print("ELPD mod.h")
    print(loo_h)
    print("ELPD mod.i")
    print(loo_i)
    print("ELPD mod.j")
    print(loo_j)
    
    print("ELPD comparison")
    print(loo_compare(loo_a_full, loo_b, loo_c, loo_d, loo_e, loo_f, #loo_g,
                      loo_h, loo_i, loo_j))
    # Model weights
    mod.a.latency$loo <- loo_a_full
    mod.b.latency$loo <- loo_b
    mod.c.latency$loo <- loo_c
    mod.d.latency$loo <- loo_d
    mod.e.latency$loo <- loo_e
    mod.h.latency$loo <- loo_h
    mod.i.latency$loo <- loo_i
    mod.j.latency$loo <- loo_j
    # mod.w <- model_weights(mod.a.latency, mod.b.latency, mod.c.latency,
    #                        mod.d.latency, mod.e.latency, mod.f.latency, #mod.g.latency,
    #                        weights = "stacking")
    mod.w <- model_weights(mod.a.latency, mod.b.latency, mod.c.latency,
                           mod.d.latency, mod.e.latency, mod.f.latency,
                           #mod.g.latency, # null model
                           mod.h.latency, mod.i.latency, mod.j.latency,
                           weights = "loo")
    print("Model weigths")
    
    print(sort(mod.w))
  }
}

## ---- Results: most supported model ----
{
  # Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
  # through posterior predictive checks
  fit <- mod.j.latency
  {
    # 1. Convergence diagnostics ----
    
    # Gelman-Rubin statistic (Rhat)
    summary(fit)$fixed  # Rhat values are part of the summary
    post <- as_draws_df(fit)
    
    # Traceplots for visual inspection
    library(bayesplot)
    set.seed(1234)
    mcmc_trace(post, regex_pars = "b_") # trace for fixed effects
    
    # 2. Posterior predictive checks -----
    
    # Posterior predictive checks using brms
    set.seed(1234)
    pp_check(fit)  # default overlay of observed vs simulated densities
    set.seed(1234)
    pp_check(fit, type = "stat")  # e.g., distribution of means
    set.seed(1234)
    pp_check(fit, type = "intervals")  # predictive intervals
    
    # 3. Residual checks ----
    
    # Extract residuals and fitted values
    d$resid <- residuals(fit)[, "Estimate"]
    d$fitted <- fitted(fit)[, "Estimate"]
    
    # Plot residuals vs fitted values to check for non-linearity
    library(ggplot2)
    ggplot(d, aes(x = fitted, y = resid)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Fitted values", y = "Residuals") +
      theme_minimal()
  }
  
  # Exploring outputs from the best candidate model
  {
    fit <- mod.j.latency
    
    # Extract posterior samples
    post <- as_draws_df(fit)
    post$sylvaticus <- -(post$b_hu_sp1 + post$b_hu_sp2 + post$b_hu_sp3)
    post$mu_sylvaticus <- post$b_hu_Intercept + post$sylvaticus
    # Compute posterior mean and 95% CrI
    mean_sylv <- mean(post$sylvaticus)
    CrI_sylv <- quantile(post$sylvaticus, probs = c(0.025, 0.975))
    
    mean_sylv
    CrI_sylv
    
    # d %>% filter(latency.time.adjusted > 0) %>% 
    #   group_by(sp) %>% 
    #   summarise(media = mean(latency.time.adjusted))
    # 
    # d %>% filter(latency.time.adjusted > 0) %>% 
    #   summarise(media = mean(latency.time.adjusted))
    
    # Extract beta (population-level effects)
    beta_cols <- grep("^b_", colnames(post), value = TRUE)
    prob_betas_below_0 <- colMeans(post[, beta_cols] < 0)
    prob_betas_above_0 <- colMeans(post[, beta_cols] > 0)
    
    # post$mu_sylvaticus <- post$b_Intercept - post$b_sp1 - post$b_sp2 - post$b_sp3
    # post_beta_sylv_above_0 <- mean(post$mu_sylvaticus > 0)
    # post_beta_sylv_below_0 <- mean(post$mu_sylvaticus < 0)
    
    post_beta_sylv_above_0 <- mean(post$sylvaticus > 0)
    post_beta_sylv_below_0 <- mean(post$sylvaticus < 0)
  }
  
  # print results
  {
    print(paste("posterior mean for sylvaticus = ", mean_sylv, sep = ""))
    print(paste("CrI for sylvaticus = ", CrI_sylv, sep = ""))
    print("prob_betas_below_0")
    print(prob_betas_below_0)
    print("prob_betas_above_0")
    print(prob_betas_above_0)
    print("post_beta_sylv_above_0")
    post_beta_sylv_above_0
    print("post_beta_sylv_below_0")
    post_beta_sylv_below_0
  }
} 

{
  # elevation * shannon interaction
  conditional_effects(fit, "elevation.std:shannon.std")
  conditional_effects(fit, "shannon.std:elevation.std")
}



# check from here, eventually

## ---- Results: second best model ----
{
  # Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
  # through posterior predictive checks
  rm(fit)
  fit <- mod.c.latency
  {
    # 1. Convergence diagnostics ----
    
    # Gelman-Rubin statistic (Rhat)
    summary(fit)$fixed  # Rhat values are part of the summary
    post <- as_draws_df(fit)
    
    # Traceplots for visual inspection
    library(bayesplot)
    set.seed(1234)
    mcmc_trace(post, regex_pars = "b_") # trace for fixed effects
    
    # 2. Posterior predictive checks -----
    
    # Posterior predictive checks using brms
    set.seed(1234)
    pp_check(fit)  # default overlay of observed vs simulated densities
    set.seed(1234)
    pp_check(fit, type = "stat")  # e.g., distribution of means
    set.seed(1234)
    pp_check(fit, type = "intervals")  # predictive intervals
    
    # 3. Residual checks ----
    
    # Extract residuals and fitted values
    d$resid <- residuals(fit)[, "Estimate"]
    d$fitted <- fitted(fit)[, "Estimate"]
    
    # Plot residuals vs fitted values to check for non-linearity
    library(ggplot2)
    ggplot(d, aes(x = fitted, y = resid)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Fitted values", y = "Residuals") +
      theme_minimal()
  }
  
  # Exploring outputs from the full model
  {
    
    # Extract posterior samples
    post <- as_draws_df(fit)
    post$sylvaticus <- -(post$b_sp1 + post$b_sp2 + post$b_sp3)
    # Compute posterior mean and 95% CrI
    mean_sylv <- mean(post$sylvaticus)
    CrI_sylv <- quantile(post$sylvaticus, probs = c(0.025, 0.975))
    
    mean_sylv
    CrI_sylv
    
    # Extract beta (population-level effects)
    beta_cols <- grep("^b_", colnames(post), value = TRUE)
    prob_betas_below_0 <- colMeans(post[, beta_cols] < 0)
    prob_betas_above_0 <- colMeans(post[, beta_cols] > 0)
    
    beta_int_elev_sp3 <- post$`b_elevation.std:sp3`
    prop_int_elev_sp3_below_0 <- mean(beta_int_elev_sp3 < 0)
    
    # # post$mu_sylvaticus <- post$b_Intercept - post$b_sp1 - post$b_sp2 - post$b_sp3
    # post_beta_sylv_above_0 <- mean(post$sylvaticus > 0)
    # post_beta_sylv_below_0 <- mean(post$sylvaticus < 0)
  }
  
  # print results
  {
    print(paste("posterior mean for sylvaticus = ", mean_sylv, sep = ""))
    print(paste("CrI for sylvaticus = ", CrI_sylv, sep = ""))
    print("prob_betas_below_0")
    print(prob_betas_below_0)
    print("prob_betas_above_0")
    print(prob_betas_above_0)
    # print("post_beta_sylv_above_0")
    # post_beta_sylv_above_0
    prop_int_elev_sp3_below_0
  }
}  