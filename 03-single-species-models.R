## ---- Single species models ----

set.seed(1234)

## ---- Loading packages ----
# devtools::install_github('bbc/bbplot')
{
  require(ggplot2)
  require(ggeffects)
  require(reshape)
  require(reshape2)
  require(lubridate)
  require(sp)
  require(car)
  # require(rgdal)
  require(plyr)
  require(dplyr)
  require(ggmap)
  register_google(key ="AIzaSyDV9QP7CdmQ3YApschivng59DoRY3rs8Rg") #È necessario attivare una chiave
  require(brms)
  require(gstat)
  require(sp)
  require(raster)
  require(foreign)
  require(car)
  require(raster)
  require(sf)
  require(ggfortify)
  require(NbClust)
  require(cluster)
  require(factoextra)
  require(ggeffects)
  require(sjPlot)
  require(performance)
  require(see)
  require(mgcv)
  require(ggeffects)
  require(dagitty)
  require(ggdag)
  require(lavaan)
  require(lmtest)
  require(ROCR)
  require(gridExtra)
  require(cluster)
  require(factoextra)
  require(psy)
  require(psych)
  require(varhandle)
  require(clValid)
  require(readxl)
  require(foreign)
  require(tibble)
  require(purrr)
  require(reshape)
  require(reshape2)
  require(lubridate)
  require(sp)
  require(sf)
  require(gstat)
  require(raster)
  require(geosphere)
  require(brms)
  require(mgcv)
  require(gam)
  require(modelr)
  require(gridExtra)
  require(randomForest)
  require(randomForestSRC)
  require(ggRandomForests)
  require(ggplotify)
  # require(ggSpatial)
  library(ncdf4) # package for netcdf manipulation
  require(plotly)
  require(bayestestR)
  require(bayesplot)
  require(fitdistrplus)
  require(hrbrthemes)
  # require(bbplot)
  require(devtools)
  require(bbplot)
  
}

## ---- Loading data ----
{
  {
    df_temp <- read_excel("data/Behaviour_habitat.xlsx") # read.csv("data/Behaviour_habitat.csv")
    head(df_temp)
  }
  
  # missing variables in other datasets
  {
    df_alp <- read_excel("data/alp.xlsx") %>% # Apodemus alpicola
      dplyr::select(ID, shannon, veg_cov)
    head(df_alp)
    
    df_fla <- read_excel("data/fla.xlsx") %>% # Apodemus flavicollis
      dplyr::select(ID, shannon, veg_cov)
    head(df_fla)
    
    df_myo <- read_excel("data/myo.xlsx") %>% # Myodes glareolus
      dplyr::select(ID, shannon, veg_cov)
    head(df_myo)
    
    df_syl <- read_excel("data/syl.xlsx") %>% # Apodemus sylvaticus
      dplyr::select(ID, shannon, veg_cov)
    head(df_syl)
    
    df_species <- bind_rows(df_alp, df_fla, df_myo, df_syl)
  }
  
  {
    # join
    df <- left_join(df_temp, df_species, join_by(id == ID))
    df <- rename(df, veg.cov = veg_cov)
  }
  
  # clean environment
  rm(df_alp, df_fla, df_myo, df_syl, df_species, df_temp)
  
  #Adjusting the date
  {
    df$date <- dmy(gsub("/", "-", df$date))
  }
  
  #Sex
  {
    df$sex <- factor(df$sex, labels=c("M", "F"))
  }
  
  
  #Selecting the variables of interest form the original dataset
  {
    d <- dplyr::select(df, date:area, sex,
                       distance.traveled, latency.time,
                       veg.cov, shannon,
                       tenv.min,
                       foot.length, tbnight:elevation) 
  }
  
  #Removing one A. sylvaticus who traveled an anomalous amount of distance
  {
    # d <- d[-which(d$id%in%d[which(d$sp=="sylvaticus" & d$distance.traveled>20000),"id"]),]
    d <- filter(d, distance.traveled < 20000)
  }
  
  # Correct Clethrionomys glareolus name
  {
    d$gen[d$gen == "Myodes"] <- "Clethrionomys"
  }
}  

#Creating four different datasets
{
  aa <- filter(d, sp=="alpicola")
  af <- filter(d, sp=="flavicollis")
  as <- filter(d, sp=="sylvaticus")
  cg <- filter(d, sp=="glareolus")
}



## ---- Exploratory analyses of response variables -----

d.sp <- list(aa, af, as, cg)

for(i in 1:length(d.sp)) {
  
  ## Distance travelled
  {
    boxplot(d.sp[[i]]$distance.traveled)
    
    #Exploring the distribution of the response variable
    {
      print(ggplot(data=d.sp[[i]]) + 
        geom_histogram(aes(x=distance.traveled), bins = 9, colour = "white", fill = "#1380A1") +
        geom_hline(yintercept = 0, size = 1, colour="#333333") +
        bbc_style() +
        labs(x="cm", y="Frequency", title=paste(d.sp[[i]]$gen[1], d.sp[[i]]$sp[1],
                                                "(n = ", dim(d.sp[[i]])[1],")", sep = " "), 
             subtitle ="Distance covered by individuals, once put in the cage") +
        scale_x_continuous(breaks=c(0, 5000, 10000, 15000), labels=c("0", "5000", "10000", "15000 cm")) +
        theme(plot.subtitle = element_text(size=18),
              plot.title = element_text(face="bold.italic")))
    }
    
    # distance travelled by the different species
    d.sp[[i]] %>% 
      summarise(mean_dist = mean(distance.traveled), sd_dist = sd(distance.traveled))
  }  
  
  ## Latency time
  {
    #Exploring the distribution of the response variables (only animals that jumped)
    {
      d.sp[[i]]$latency.time.adjusted <- d.sp[[i]]$latency.time
      d.sp[[i]][which(d.sp[[i]]$latency.time>=570), "latency.time.adjusted"] <- 0
      hist(d.sp[[i]]$latency.time.adjusted)
      print(ggplot(data=d.sp[[i]][which(d.sp[[i]]$latency.time.adjusted!=0),]) + 
        geom_histogram(aes(x=latency.time.adjusted), bins = 9, colour = "white", fill = "#1380A1") +
        geom_hline(yintercept = 0, size = 1, colour="#333333") +
        bbc_style() +
        labs(x="seconds", y="Frequency", title=paste(d.sp[[i]]$gen[1], d.sp[[i]]$sp[1],
                                                     "(n = ", dim(d.sp[[i]])[1],")", sep = " "),
             subtitle="Time before individuals jumped, once put in the cage (non-zero)") +
        scale_x_continuous(breaks=c(0, 200, 400, 600), labels=c("0", "200", "400", "600 s")) +
        theme(plot.subtitle = element_text(size=18)))
    }
  }
}

# Define priors
{
  ## priors for single species models, with random effect
  # priors.rf <- c(
  #   # Priors for the intercept
  #   prior(normal(7, 2), class = "Intercept"),
  #   # Prior for the shape parameter of the Gamma distribution
  #   prior(gamma(0.01, 0.01), class = "shape"),  # or exponential(1)
  #   # Prior for the group-level standard deviation (random effect of area)
  #   prior(exponential(1), class = "sd", group = "area")
  # )

  # priors for single species models, without random effect
  priors.no.rf <- c(
    # Priors for the intercept
    prior(normal(7, 2), class = "Intercept"),
    # Prior for the shape parameter of the Gamma distribution
    prior(gamma(0.01, 0.01), class = "shape")  # or exponential(1)
  )
}

dataset.sp <- list(aa, af, as, cg)
var_draws <- list()



## ---- Models for single species - distance travelled ----
{
  # Fit models and extract poterior for variance
  for (i in 1:length(dataset.sp)) {
    
    # Fit models
    {
      # mod.distance.rf <- brm(distance.traveled ~ 1 + (1 | area),
      #                        family=Gamma(link="log"), data=dataset,
      #                        chains=4, cores=4, iter=5000, warmup=2000, #thin = 10,
      #                        control = list(adapt_delta = 0.9999, max_treedepth = 15),
      #                        prior = priors.rf,
      #                        seed = 1234)
      mod.distance.no.rf <- brm(distance.traveled ~ 1,
                                family=Gamma(link="log"),
                                data=dataset.sp[[i]],
                                chains=4, cores=4, iter=5000, warmup=2000, #thin = 10,
                                control = list(adapt_delta = 0.9999, max_treedepth = 15),
                                prior = priors.no.rf,
                                seed = 1234)
    }
    
    #Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
    #through posterior predictive checks
    {
      # summary(mod.distance.rf)
      # brms::pp_check(mod.distance.rf)
      # plot(mod.distance.rf)
      
      summary(mod.distance.no.rf)
      brms::pp_check(mod.distance.no.rf) # if the two pp_checks look similar, the variance absorbed by area is not relevan
      plot(mod.distance.no.rf)
    }
    
    # Posterior distribution of variance
    {
      post <- as_draws_df(mod.distance.no.rf)
      # Get posterior draws of the intercept and shape
      mu_draws <- exp(post$b_Intercept)
      shape_draws <- post$shape
      # Compute posterior distribution of variance
      var_draws[[i]] <- (mu_draws^2) / shape_draws
    }      
  }      
}

# Comparison of variances across species
{
  df_var <- data.frame(
    glareolus = var_draws[[4]],
    flavicollis = var_draws[[2]],
    sylvaticus = var_draws[[3]],
    alpicola = var_draws[[1]]
  )
  
  apply(df_var, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  # pairwise comparisons
  mean(df_var$glareolus > df_var$alpicola)  # e.g., P(glareolus variance > alpicola variance)
  
  library(ggplot2)
  library(tidyr)
  
  df_long <- pivot_longer(df_var, everything(), names_to = "species", values_to = "variance")
  df_long_distance <- df_long
  
  species_colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
  
  ggplot(df_long, aes(x = variance, fill = species)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = species_colors) +
    theme_minimal() +
    xlab("Variance in distance travelled") +
    ylab("Density") +
    labs(title = "") -> fig4a
  print(fig4a)
  
  # Example: probability that variance of species A > species B
  mean(df_var$sylvaticus > df_var$alpicola)       # P(sylvaticus > alpicola)
  mean(df_var$sylvaticus > df_var$glareolus)      # etc.
  mean(df_var$sylvaticus > df_var$flavicollis) 
  mean(df_var$glareolus > df_var$alpicola)
  mean(df_var$glareolus > df_var$flavicollis)
  mean(df_var$flavicollis > df_var$alpicola)
  
  species <- colnames(df_var)
  n <- length(species)
  
  prob_matrix <- matrix(NA, n, n, dimnames = list(species, species))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        prob_matrix[i, j] <- mean(df_var[[species[i]]] > df_var[[species[j]]])
      }
    }
  }
  
  round(prob_matrix, 3)
  
}


## ---- Models for single species - latency time ----
{
  rm(mod.distance.no.rf)
  # Fit models and extract poterior for variance
  for (i in 1:length(dataset.sp)) {
    
    # Fit models
    {
      # mod.distance.rf <- brm(distance.traveled ~ 1 + (1 | area),
      #                        family=Gamma(link="log"), data=dataset,
      #                        chains=4, cores=4, iter=5000, warmup=2000, #thin = 10,
      #                        control = list(adapt_delta = 0.9999, max_treedepth = 15),
      #                        prior = priors.rf,
      #                        seed = 1234)
      mod.distance.no.rf <- brm(latency.time ~ 1,
                                family=Gamma(link="log"),
                                
                                ## NB, change the dataset to consider only jumpers:
                                # data= (dataset.sp[[i]])[which(dataset.sp[[i]]$latency.time!=570),],
                                
                                ## or all individuals (censored variable):
                                data = dataset.sp[[i]],
                                
                                chains=4, cores=4, iter=5000, warmup=2000, #thin = 10,
                                control = list(adapt_delta = 0.9999, max_treedepth = 15),
                                prior = priors.no.rf,
                                seed = 1234)
    }
    
    #Exploring convergence, the posterior distribution of model parameters and goodness-of-fit
    #through posterior predictive checks
    {
      # summary(mod.distance.rf)
      # brms::pp_check(mod.distance.rf)
      # plot(mod.distance.rf)
      
      summary(mod.distance.no.rf)
      brms::pp_check(mod.distance.no.rf) # if the two pp_checks look similar, the variance absorbed by area is not relevan
      plot(mod.distance.no.rf)
    }
    
    # Posterior distribution of variance
    {
      post <- as_draws_df(mod.distance.no.rf)
      # Get posterior draws of the intercept and shape
      mu_draws <- exp(post$b_Intercept)
      shape_draws <- post$shape
      # Compute posterior distribution of variance
      var_draws[[i]] <- (mu_draws^2) / shape_draws
    }      
  }      
}

# Comparison of variances across species
{
  df_var <- data.frame(
    glareolus = var_draws[[4]],
    flavicollis = var_draws[[2]],
    sylvaticus = var_draws[[3]],
    alpicola = var_draws[[1]]
  )
  
  apply(df_var, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
  # pairwise comparisons
  mean(df_var$glareolus > df_var$alpicola)  # e.g., P(glareolus variance > alpicola variance)
  
  library(ggplot2)
  library(tidyr)
  
  df_long <- pivot_longer(df_var, everything(), names_to = "species", values_to = "variance")
  df_long_latency <- df_long
  
  ggplot(df_long, aes(x = variance, fill = species)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = species_colors) +
    theme_minimal() +
    xlab("Variance in latency time") +
    ylab("Density") +
    labs(title = "") -> fig4b
  print(fig4b)
  
  # Example: probability that variance of species A > species B
  mean(df_var$sylvaticus > df_var$alpicola)       # P(sylvaticus > alpicola)
  mean(df_var$sylvaticus > df_var$glareolus)      # etc.
  mean(df_var$sylvaticus > df_var$flavicollis) 
  mean(df_var$glareolus > df_var$alpicola)
  mean(df_var$glareolus > df_var$flavicollis)
  mean(df_var$flavicollis > df_var$alpicola)
  
  species <- colnames(df_var)
  n <- length(species)
  
  prob_matrix <- matrix(NA, n, n, dimnames = list(species, species))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        prob_matrix[i, j] <- mean(df_var[[species[i]]] > df_var[[species[j]]])
      }
    }
  }
  
  round(prob_matrix, 3)
  
}

## ---- plot for paper fig4 ----

# Create a mapping from original names to cleaned names
{
  species_mapping <- c(
    "alpicola" = "A. alpicola",
    "flavicollis" = "A. flavicollis",
    "glareolus" = "C. glareolus",
    "sylvaticus"   = "A. sylvaticus"
  )
  
  # Add a new column to your data with the cleaned species names
  df_long_distance$Species <- species_mapping[df_long_distance$species]
  df_long_latency$Species <- species_mapping[df_long_latency$species]
  
  # save(df_long_distance, file = "output/df_long_distance.RData", compress = F)
  # save(df_long_latency, file = "output/df_long_latency.RData", compress = F)
}

# Use italics for species
{
  italic_labels <- setNames(
    lapply(df_long_distance$Species, function(sp) bquote(italic(.(gsub("_", " ", sp))))),
    df_long_distance$Species
  )
}

## ---- Colors for plots ----
{
  species_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
}

space <- scales::format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

## Plots
{
  # Variance in distance travelled
  {
    ggplot(df_long_distance, aes(x = variance, fill = Species, col = Species)) +
      geom_density(alpha = 0.6, col = NA) +
      scale_color_manual(values = species_colors) +
      scale_fill_manual(
        values = species_colors,
        labels = italic_labels
      ) +
      theme_minimal() +
      xlab("Variance in distance travelled") +
      ylab("Density") +
      # scale_x_continuous(labels = scales::comma_format()) +
      scale_x_continuous(labels = space, limits = c(0,40000000)) +
      # scale_y_continuous(labels = space) +
      guides(color = "none") +
      labs(title = "") -> fig4a
    print(fig4a)
  }
  
  # Variance in latency time
  {
    ggplot(df_long_latency, aes(x = variance, fill = Species, col = Species)) +
      geom_density(alpha = 0.6, col = NA) +
      scale_color_manual(values = species_colors) +
      scale_fill_manual(
        values = species_colors,
        labels = italic_labels
      ) +
      theme_minimal() +
      xlab("Variance in latency time") +
      ylab("Density") +
      # scale_x_continuous(labels = scales::comma_format()) +
      scale_x_continuous(labels = space, limits = c(0,400000)) +
      # scale_y_continuous(labels = space) +
      guides(color = "none") +
      labs(title = "") -> fig4b
    print(fig4b)
  }
  
  # Combined plot
  {
    library(patchwork)  
    (fig4a / fig4b) + 
      plot_layout(guides = "collect") +          # <-- collect legends here
      plot_annotation(tag_levels = 'a') & 
      theme(legend.position = "bottom") -> fig4
    print(fig4)
  }
  
  # Save combined plot
  {
    # tiff("figs/fig4.tiff", res=1000, width = 18, height = 14, units = "cm")
    # print(fig4)
    # dev.off()
  }
}

# Alternative Ridge Plots
{
  # Distance travelled
  {
    ggplot(df_long_distance, aes(x = variance, y = Species, fill = Species)) +
      geom_density_ridges(scale = 1.2, alpha = 0.6, rel_min_height = 0.01, col = NA) +
      # scale_color_manual(values = species_colors) +
      scale_fill_manual(
        values = species_colors,
        labels = italic_labels
      ) +
      theme_minimal() +
      xlab("Variance in distance travelled") +
      ylab("Density") +
      # scale_x_continuous(labels = scales::comma_format()) +
      scale_x_continuous(labels = space, limits = c(0,30000000)) +
      scale_y_discrete(labels = italic_labels) +
      theme(legend.position = "none") +
      labs(title = "") -> fig4a
    print(fig4a)
  }
  
  # Latency time
  {
    ggplot(df_long_latency, aes(x = variance, y = Species, fill = Species)) +
      geom_density_ridges(scale = 1.2, alpha = 0.6, rel_min_height = 0.01, col = NA) +
      # scale_color_manual(values = species_colors) +
      scale_fill_manual(
        values = species_colors,
        labels = italic_labels
      ) +
      theme_minimal() +
      xlab("Variance in latency time") +
      ylab("") +
      # scale_x_continuous(labels = scales::comma_format()) +
      scale_x_continuous(labels = space, limits = c(0,300000)) +
      scale_y_discrete(labels = NULL) +
      theme(legend.position = "none") +
      labs(title = "") -> fig4b
    print(fig4b)
  }

  # Combined plot
  {
    library(patchwork)  
    (fig4a + fig4b) + 
      plot_layout(guides = "collect") +          # <-- collect legends here
      plot_annotation(tag_levels = 'a') & 
      theme(legend.position = "bottom") -> fig4_v2
    print(fig4_v2)
  }
  
  # Save combined plot
  {
    # tiff("figs/fig4_v2.tiff", res=1000, width = 18, height = 14, units = "cm")
    # print(fig4_v2)
    # dev.off()
  }
}