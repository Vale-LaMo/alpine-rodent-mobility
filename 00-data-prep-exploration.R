# Effect of environmental and individual factors

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
  library(tidyr)
  
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
}  
  
## ---- Exploratory analyses of response variables -----

## Distance travelled
{
  boxplot(d$distance.traveled)
  
  #Removing one A. sylvaticus who traveled an anomalous amount of distance
  {
    # d <- d[-which(d$id%in%d[which(d$sp=="sylvaticus" & d$distance.traveled>20000),"id"]),]
    d <- filter(d, distance.traveled < 20000)
  }
  
  #Exploring the distribution of the response variable
  {
    ggplot(data=d) + 
      geom_histogram(aes(x=distance.traveled), bins = 9, colour = "white", fill = "#1380A1") +
      geom_hline(yintercept = 0, size = 1, colour="#333333") +
      bbc_style() +
      labs(x="cm", y="Frequency", title=paste("All species (n = ", dim(d)[1],")", sep = ""), 
           subtitle ="Distance covered by individuals, once put in the cage") +
      scale_x_continuous(breaks=c(0, 5000, 10000, 15000), labels=c("0", "5000", "10000", "15000 cm")) +
      theme(plot.subtitle = element_text(size=18),
            plot.title = element_text(face="bold.italic"))
  }
  
  # distance travelled by the different species
  d %>% group_by(sp) %>% 
    summarise(mean_dist = mean(distance.traveled), sd_dist = sd(distance.traveled))
}  

## Latency time
{
  #Exploring the distribution of the response variables (only animals that jumped)
  {
    d$latency.time.adjusted <- d$latency.time
    d[which(d$latency.time>=570), "latency.time.adjusted"] <- 0
    hist(d$latency.time.adjusted)
    ggplot(data=d[which(d$latency.time.adjusted!=0),]) + 
      geom_histogram(aes(x=latency.time.adjusted), bins = 9, colour = "white", fill = "#1380A1") +
      geom_hline(yintercept = 0, size = 1, colour="#333333") +
      bbc_style() +
      labs(x="seconds", y="Frequency", title=paste("All species (n = ", dim(d)[1],")", sep = ""),
           subtitle="Time before individuals jumped, once put in the cage (non-zero)") +
      scale_x_continuous(breaks=c(0, 200, 400, 600), labels=c("0", "200", "400", "600 s")) +
      theme(plot.subtitle = element_text(size=18))
  }
  
  #Comparing the proportion of generated zeroes between species in latency time
  {
    d[which(d$latency.time.adjusted==0), "latency.zero"] <- 0
    d[which(d$latency.time.adjusted!=0), "latency.zero"] <- 1
    round(prop.table(xtabs(~latency.zero + sp, data=d), margin=2), 3)
  }
    
  # fig3
  {
    # Step 1: Create the table of proportions
    tab <- xtabs(~latency.zero + sp, data = d)
    tab_df <- as.data.frame.matrix(tab)
    colnames(tab_df) <- c("A. alpicola", "A. flavicollis", "C. glareolus", "A. sylvaticus")
    
    # Step 2: Convert to long format
    tab_df$latency.zero <- rownames(tab_df)
    tab_long <- pivot_longer(tab_df, cols = -latency.zero, names_to = "species", values_to = "no. individuals")
    
    # Step 3: Make labels more informative (optional)
    tab_long$latency.zero <- factor(tab_long$latency.zero, 
                                    levels = c("0", "1"),
                                    labels = c("Did Not Jump", "Jumped"))
    
    # Step 4: Plot
    library(ggokabeito)
    
    # Use italics for species
    {
      italic_labels <- setNames(
        lapply(tab_long$species, function(sp) bquote(italic(.(gsub("_", " ", sp))))),
        tab_long$species
      )
    }
    
    library(see)
    see::oi_colors()
    # #    orange light blue      green      amber       blue        red     purple
    # "#E69F00"  "#56B4E9"  "#009E73"  "#F5C710"  "#0072B2"  "#D55E00"  "#CC79A7"
    # grey      black
    # "#999999"  "#000000"
    
    # save(tab_long, file = "output/tab_long.RData", compress = F)
    ggplot(tab_long, aes(x = species, y = `no. individuals`, fill = latency.zero)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = c("#F5C710", "#0072B2")) +
      scale_x_discrete(labels = italic_labels) +
      # can be one of "full", "black_first", full_original, or black_original. 
      labs(
        x = "Species",
        y = "No. individuals",
        fill = "Escape behaviour",
        title = ""
      ) +
      theme_minimal(base_size = 14) -> fig3
    fig3
    
    # tiff("figs/fig3.tiff", res=1000, width = 18, height = 14, units = "cm")
    # print(fig3)
    # dev.off()
  }
}


## ---- Correlation among behaviours ----

# Summary statics about behaviours
{
  summary(d$distance.traveled) # range distance travelled (excluded one oulier reported above)
  d %>% filter(latency.time < 570) %>% 
    summarise(min=min(latency.time), max=max(latency.time)) # range latency time for animals that jumped
  d %>% filter(latency.time >= 570) %>% 
    count(latency.time) # number of animal that did not jump
  d %>% filter(latency.time < 570) %>% 
    dim() # number of animal that did jump
  # 142+59 check
  plot(d$distance.traveled, d$latency.time)
}

# Analyses of correlation using all animals
{
  # Correlation model
  {
    # Rescale if needed to improve model convergence
    d$distance_scaled <- scale(d$distance.traveled)
    d$latency_scaled <- scale(d$latency.time)
    
    # Multivariate model: jointly model both variables with correlation
    model_corr <- brm(
      mvbind(distance_scaled, latency_scaled) ~ 1,  # no predictors, just correlation
      data = d,
      family = gaussian(),
      # rescor = TRUE,
      chains = 4,
      cores = 4,
      iter = 5000,
      warmup = 2000,
      save_pars = save_pars(all = TRUE),
      backend = "cmdstanr",
      seed = 1234
    ) #+
    # set_rescor(rescor = TRUE)
  }
  
  # View summary
  {
    print(summary(model_corr))
    # If the 95% Credible Interval does not overlap 0, then you have good evidence for a correlation.
    draws <- as_draws_df(model_corr) # extract posteriors
    colnames(draws)
    cor_samples <- draws$rescor__distancescaled__latencyscaled
    prob_corr_lt_0 <- mean(cor_samples < 0) # Calculate the posterior probability that the correlation < 0
    print(paste("posterior probability that the correlation < 0 = ", prob_corr_lt_0, sep = ""))
    # Also get the posterior mean and credible interval
    mean(cor_samples)
    quantile(cor_samples, probs = c(0.025, 0.975)) # coincide con output del summary
  }
}

# Analyses of correlation using only jumpers    
{
  # Correlation model
  {
    d %>% filter(d$latency.time < 570) -> d.jumpers
    
    # Rescale if needed to improve model convergence
    d.jumpers$distance_scaled <- scale(d.jumpers$distance.traveled)
    d.jumpers$latency_scaled <- scale(d.jumpers$latency.time)
    
    # Multivariate model: jointly model both variables with correlation
    model_corr <- brm(
      mvbind(distance_scaled, latency_scaled) ~ 1,  # no predictors, just correlation
      data = d.jumpers,
      family = gaussian(),
      # rescor = TRUE,
      chains = 4,
      cores = 4,
      iter = 5000,
      warmup = 2000,
      save_pars = save_pars(all = TRUE),
      backend = "cmdstanr",
      seed = 1234
    ) #+
    # set_rescor(rescor = TRUE)
  }
  
  # View summary
  {
    print(summary(model_corr))
    # If the 95% Credible Interval does not overlap 0, then you have good evidence for a correlation.
    draws <- as_draws_df(model_corr) # extract posteriors
    colnames(draws)
    cor_samples <- draws$rescor__distancescaled__latencyscaled
    prob_corr_lt_0 <- mean(cor_samples < 0) # Calculate the posterior probability that the correlation < 0
    print(paste("posterior probability that the correlation < 0 = ", prob_corr_lt_0, sep = ""))
    # Also get the posterior mean and credible interval
    mean(cor_samples)
    quantile(cor_samples, probs = c(0.025, 0.975)) # coincide con output del summary
  }  
}

## ---- Exploratory analyses of predictors -----

## Altitudinal distribution of species
{
  boxplot(d$elevation ~ d$sp)
  ggplot(d) + geom_histogram(aes(elevation)) + facet_grid(vars(sp))
  
  cut(d$elevation, breaks = 6) # intervals of about 250 m
  d$elevation.int <- cut(d$elevation, breaks = 6) 
  table(d$sp, d$elevation.int)
  colSums(table(d$sp, d$elevation.int))
  
  d %>% 
    group_by(sp) %>% 
    summarise(min = min(elevation), max = max(elevation))
}
  

## ---- Set-up for models ----
{
  # Standardization of predictors
  d$elevation.std <- (d$elevation - mean(d$elevation))/sd(d$elevation)
  d$shannon.std <- (d$shannon - mean(d$shannon))/sd(d$shannon)
  
  #Setting sum-to-zero contrasts for species
  {
    d$sp <- factor(d$sp)
    contrasts(d$sp) <- contr.sum
  }
}

# save(d, file = "output/d.RData", compress = FALSE)
