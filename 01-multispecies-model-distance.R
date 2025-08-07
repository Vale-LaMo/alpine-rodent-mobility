# Modelling species-specific differences on traveled distance
library(tidyr)

set.seed(1234)
source("00-data-prep-exploration.R")

## ---- Model fitting ----
{
  # Gamma with species interactions, full model
  
  # Define priors
  {
    priors <- c(
      # Priors for the intercept
      prior(normal(7, 2), class = "Intercept"),
      # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
      prior(normal(0, 2), class = "b"),
      # Prior for the shape parameter of the Gamma distribution
      prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
      # Prior for the group-level standard deviation (random effect of area)
      prior(exponential(1), class = "sd", group = "area")
    )
  }
  
  # Fit full model
  {
    mod.a.distance <- brm(distance.traveled ~ elevation.std*sp + shannon.std*sp + sex +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=priors,
                          seed = 1234)
    # save(mod.a.distance, file = "output/mod.a.distance.RData", compress = FALSE)
  }
  
  # Multicollinearity check
  {
    library(posterior)
    # Extract posterior draws as a data frame
    post <- as_draws_df(mod.a.distance)
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
      # Priors for the intercept
      prior(normal(7, 2), class = "Intercept"),
      # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
      prior(normal(0, 5), class = "b"),
      # Prior for the shape parameter of the Gamma distribution
      prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
      # Prior for the group-level standard deviation (random effect of area)
      prior(exponential(1), class = "sd", group = "area")
    )
    mod.a.distance.sens <- brm(distance.traveled ~ elevation.std*sp + shannon.std*sp + sex +
                                 (1 | area),
                               family=Gamma(link="log"),
                               data=d,
                               chains=4, cores=4, 
                               iter=5000, warmup=2000,
                               # thin = 10,
                               control = list(adapt_delta = 0.9999, max_treedepth = 15),
                               prior=priors.sens)
    # save(mod.a.distance.sens, file = "output/mod.a.distance.sens.RData", compress = FALSE)
  }
  
  # print results
  {
    print(summary(mod.a.distance))
    print(summary(mod.a.distance.sens))
    print(posterior_cor[elevation_terms, shannon_terms])
  }
}


## ---- Model selection ----
{
  # Comparing models through leave-one-out cross validation
  
  # Fitting nested models
  {
    mod.b.distance <- brm(distance.traveled ~ elevation.std*sp +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=priors,
                          seed = 1234)
    mod.c.distance <- brm(distance.traveled ~ elevation.std + sp +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=priors,
                          seed = 1234)
    mod.d.distance <- brm(distance.traveled ~ sp +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=priors,
                          seed = 1234)
    mod.e.distance <- brm(distance.traveled ~ 1 +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=c(
                            # Priors for the intercept
                            prior(normal(7, 2), class = "Intercept"),
                            # # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
                            # prior(normal(0, 2), class = "b"),
                            # Prior for the shape parameter of the Gamma distribution
                            prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
                            # Prior for the group-level standard deviation (random effect of area)
                            prior(exponential(1), class = "sd", group = "area")),
                          seed = 1234)
    # mod.f.distance <- brm(distance.traveled ~ shannon.std*sp +
    #                         (1 | area),
    #                       family=Gamma(link="log"),
    #                       data=d,
    #                       chains=4, cores=4, 
    #                       iter=5000, warmup=2000,
    #                       # thin = 10,
    #                       control = list(adapt_delta = 0.9999, max_treedepth = 15),
    #                       prior=priors,
    #                       seed = 1234)
    # mod.g.distance <- brm(distance.traveled ~ shannon.std + sp +
    #                         (1 | area),
    #                       family=Gamma(link="log"),
    #                       data=d,
    #                       chains=4, cores=4, 
    #                       iter=5000, warmup=2000,
    #                       # thin = 10,
    #                       control = list(adapt_delta = 0.9999, max_treedepth = 15),
    #                       prior=priors
    #                       seed = 1234)
    mod.h.distance <- brm(distance.traveled ~ elevation.std * shannon.std + sp + sex +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=c(
                            # Priors for the intercept
                            prior(normal(7, 2), class = "Intercept"),
                            # # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
                            prior(normal(0, 2), class = "b"),
                            # Prior for the shape parameter of the Gamma distribution
                            prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
                            # Prior for the group-level standard deviation (random effect of area)
                            prior(exponential(1), class = "sd", group = "area")),
                          seed = 1234)
    mod.i.distance <- brm(distance.traveled ~ elevation.std * shannon.std * sp + sex +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=c(
                            # Priors for the intercept
                            prior(normal(7, 2), class = "Intercept"),
                            # # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
                            prior(normal(0, 2), class = "b"),
                            # Prior for the shape parameter of the Gamma distribution
                            prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
                            # Prior for the group-level standard deviation (random effect of area)
                            prior(exponential(1), class = "sd", group = "area")),
                          seed = 1234)
    mod.j.distance <- brm(distance.traveled ~ elevation.std * shannon.std +
                            (1 | area),
                          family=Gamma(link="log"),
                          data=d,
                          chains=4, cores=4, 
                          iter=5000, warmup=2000,
                          # thin = 10,
                          control = list(adapt_delta = 0.9999, max_treedepth = 15),
                          prior=c(
                            # Priors for the intercept
                            prior(normal(7, 2), class = "Intercept"),
                            # # Priors for slopes (standardized predictors => ~N(0, 1) or N(0, 2) is common)
                            prior(normal(0, 2), class = "b"),
                            # Prior for the shape parameter of the Gamma distribution
                            prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
                            # Prior for the group-level standard deviation (random effect of area)
                            prior(exponential(1), class = "sd", group = "area")),
                          seed = 1234)
    
    # save(mod.b.distance, file = "output/mod.b.distance.RData", compress = FALSE)
    # save(mod.c.distance, file = "output/mod.c.distance.RData", compress = FALSE)
    # save(mod.d.distance, file = "output/mod.d.distance.RData", compress = FALSE)
    # save(mod.e.distance, file = "output/mod.e.distance.RData", compress = FALSE)
    # save(mod.h.distance, file = "output/mod.h.distance.RData", compress = FALSE)
    # save(mod.i.distance, file = "output/mod.i.distance.RData", compress = FALSE)
    # save(mod.j.distance, file = "output/mod.j.distance.RData", compress = FALSE)
  }
  
  # LOO cross-validation
  # load("output/mod.a.distance.RData")
  # load("output/mod.b.distance.RData")
  # load("output/mod.c.distance.RData")
  # load("output/mod.d.distance.RData")
  # load("output/mod.e.distance.RData")
  # load("output/mod.h.distance.RData")
  # load("output/mod.i.distance.RData")
  # load("output/mod.j.distance.RData")
  {
    loo_a_full <- loo(mod.a.distance, reloo=TRUE, seed = 1234) # Use it only if k > 0.7 warnings in loo() output
    loo_b <- loo(mod.b.distance, reloo = TRUE, seed = 1234)
    loo_c <- loo(mod.c.distance, reloo = TRUE, seed = 1234)
    loo_d <- loo(mod.d.distance, seed = 1234)
    loo_e <- loo(mod.e.distance, seed = 1234)
    # loo_f <- loo(mod.f.distance, seed = 1234)
    # 
    loo_h <- loo(mod.h.distance, reloo = TRUE, seed = 1234)
    loo_i <- loo(mod.i.distance, reloo = TRUE, seed = 1234)
    loo_j <- loo(mod.j.distance, seed = 1234)
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
    print("ELPD mod.h")
    print(loo_h)
    print("ELPD mod.i")
    print(loo_i)
    print("ELPD mod.j")
    print(loo_j)
    
    print("ELPD comparison")
    # loo_compare(loo_a_full, loo_b, loo_c, loo_d, loo_e, loo_f, loo_g)
    print(loo_compare(loo_a_full, loo_b, loo_c, loo_d, loo_h, loo_i, loo_j))
    # Model weights
    mod.a.distance$loo <- loo_a_full
    mod.b.distance$loo <- loo_b
    mod.c.distance$loo <- loo_c
    mod.d.distance$loo <- loo_d
    mod.e.distance$loo <- loo_e
    mod.h.distance$loo <- loo_h
    mod.i.distance$loo <- loo_i
    mod.j.distance$loo <- loo_j
    # mod.w <- model_weights(mod.a.distance, mod.b.distance, mod.c.distance,
    #                        mod.d.distance, mod.h.distance, mod.i.distance,
    #                        weights = "stacking")
    mod.w <- model_weights(mod.a.distance, mod.b.distance, mod.c.distance,
                           mod.d.distance, mod.h.distance, mod.i.distance, mod.j.distance,
                           weights = "loo")
    print("Model weigths")
    print(sort(mod.w))
  }
}

## ---- Results: most supported model ----
{
  # Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
  # through posterior predictive checks
  fit <- mod.d.distance
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
    fit <- mod.d.distance
    
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
    
    # post$mu_sylvaticus <- post$b_Intercept - post$b_sp1 - post$b_sp2 - post$b_sp3
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
  }
  
  # fig2  
  {
    library(tidybayes)
    library(tidyr)
    
    # Compute log-scale estimates for each species
    log_preds <- post %>%
      transmute(
        'A. alpicola'= b_Intercept + b_sp1,
        'A. flavicollis' = b_Intercept + b_sp2,
        'C. glareolus'   = b_Intercept + b_sp3,
        'A. sylvaticus'  = b_Intercept - (b_sp1 + b_sp2 + b_sp3)  # due to sum-to-zero contrasts
      ) %>%
      pivot_longer(cols = everything(), names_to = "sp", values_to = "log_estimate")
    
    # Convert to original scale (distance in cm)
    log_preds <- log_preds %>%
      mutate(estimate = exp(log_estimate))
    
    # Summarise: mean and 95% CrI per species
    summary_sp <- log_preds %>%
      group_by(sp) %>%
      mean_qi(estimate, .width = 0.95)
    
    # Create a mapping from original names to cleaned names
    species_mapping <- c(
      "alpicola" = "A. alpicola",
      "flavicollis" = "A. flavicollis",
      "glareolus" = "C. glareolus",
      "sylvaticus"   = "A. sylvaticus"
    )
    
    # Add a new column to your data with the cleaned species names
    d$sp_clean <- species_mapping[d$sp]
    
    # Grand mean with CI (from intercept)
    grand_mean <- fixef(mod.d.distance)["Intercept", ]
    grand_mean_df <- data.frame(
      mean = exp(grand_mean["Estimate"]),
      lwr = exp(grand_mean["Q2.5"]),
      upr = exp(grand_mean["Q97.5"])
    )
    
    ## ---- Colors for plots ----
    {
      species_colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
    }
    # #    orange light blue      green      amber       blue        red     purple
    # "#E69F00"  "#56B4E9"  "#009E73"  "#F5C710"  "#0072B2"  "#D55E00"  "#CC79A7"
    # grey      black
    # "#999999"  "#000000"
    
    # save(summary_sp, file = "output/summary_sp.RData", compress = FALSE)
    fig2 <- ggplot(summary_sp, aes(x = sp, y = estimate, col = sp)) +
      geom_jitter(data = d, aes(x = sp_clean, y = distance.traveled),
                  width = 0.15, alpha = 0.2, color = "black") +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = grand_mean_df$lwr, ymax = grand_mean_df$upr),
                fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
      geom_point(size = 2, col = species_colors) +
      geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.1, lwd=1, col = species_colors) +
      labs(
        x = "Species",
        y = "Estimated distance travelled (cm)",
        # title = "Posterior estimates and raw data per species"
      ) +
      theme_minimal() +
      geom_hline(yintercept = grand_mean_df$mean, linetype = "dashed", color = "grey")
    fig2
    
    # tiff("figs/fig2.tiff", res=1000, width = 18, height = 14, units = "cm")
    # print(fig2)
    # dev.off()
  }
  fig2
}  

## ---- Results: second most supported model ----
{
  # Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
  # through posterior predictive checks
  rm(fit)
  fit <- mod.c.distance
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
  
  # Exploring outputs from the model
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
    
    post$mu_sylvaticus <- post$b_Intercept - post$b_sp1 - post$b_sp2 - post$b_sp3
    post_beta_sylv_above_0 <- mean(post$mu_sylvaticus > 0)
    post_beta_sylv_below_0 <- mean(post$mu_sylvaticus < 0)
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
  }
}  


## ---- Results: full model ----
{
  # Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
  # through posterior predictive checks
  rm(fit)
  fit <- mod.a.distance
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
    
    post$mu_sylvaticus <- post$b_Intercept - post$b_sp1 - post$b_sp2 - post$b_sp3
    post_beta_sylv_above_0 <- mean(post$mu_sylvaticus > 0)
    post_beta_sylv_below_0 <- mean(post$mu_sylvaticus < 0)
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
  }
}  

## ---- Results: elevation*Shannon interaction. ----
{
  conditional_effects(mod.j.distance, "elevation.std:shannon.std")
  conditional_effects(mod.j.distance, "shannon.std:elevation.std")
}