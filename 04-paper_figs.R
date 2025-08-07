## Stand-alone script for paper figs

{
  # load baseline data
  load("output/d.RData")
}

{
  # Common settings for plots
  
  # library(see)
  # see::oi_colors()
  # #    orange light blue      green      amber       blue        red     purple
  # "#E69F00"  "#56B4E9"  "#009E73"  "#F5C710"  "#0072B2"  "#D55E00"  "#CC79A7"
  # grey      black
  # "#999999"  "#000000"
  species_colors <- c("#0072B2", "#E69F00", "#009E73", "#D55E00")
  
  space <- scales::format_format(big.mark = " ", decimal.mark = ",", scientific = FALSE)

  # Create a mapping from original names to cleaned names
  {
    species_mapping <- c(
    "alpicola" = "A. alpicola",
    "flavicollis" = "A. flavicollis",
    "glareolus" = "C. glareolus",
    "sylvaticus"   = "A. sylvaticus"
  )
    
    # Add a new column to your data with the cleaned species names
    d$sp_clean <- species_mapping[d$sp]
  }
}


## ---- Fig 2 ----

{
  # load data
  load("output/summary_sp.RData")
  
  load("output/mod.d.distance.RData")
  # Grand mean with CI (from intercept)
  grand_mean <- fixef(mod.d.distance)["Intercept", ]
  grand_mean_df <- data.frame(
    mean = exp(grand_mean["Estimate"]),
    lwr = exp(grand_mean["Q2.5"]),
    upr = exp(grand_mean["Q97.5"])
  )
}

{
  # Plot settings
  # Use italics for species
  {
    italic_labels <- setNames(
      lapply(summary_sp$sp, function(sp) bquote(italic(.(gsub("_", " ", sp))))),
      summary_sp$sp
    )
  }
}


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
  scale_x_discrete(labels = italic_labels) +
  theme_minimal() +
  geom_hline(yintercept = grand_mean_df$mean, linetype = "dashed", color = "grey")
fig2
# tiff("figs/fig2.tiff", res=1000, width = 18, height = 14, units = "cm")
# print(fig2)
# dev.off()


## ---- Fig 3 ----

load("output/tab_long.RData")
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


## ---- Fig 4 ----

load("output/df_long_distance.RData")
load("output/df_long_latency.RData")

# Use italics for species
{
  italic_labels <- setNames(
    lapply(df_long_distance$Species, function(sp) bquote(italic(.(gsub("_", " ", sp))))),
    df_long_distance$Species
  )
}


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
      ggridges::geom_density_ridges(scale = 1.2, alpha = 0.6, rel_min_height = 0.01, col = NA) +
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
      ggridges::geom_density_ridges(scale = 1.2, alpha = 0.6, rel_min_height = 0.01, col = NA) +
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


##---- Fig.4 Statistics for the text ----
df_long_distance %>% 
  group_by(Species) %>% 
  summarise(median = median(variance),
            IQR = IQR(variance))

df_long_latency %>% 
  group_by(Species) %>% 
  summarise(median = median(variance),
            IQR = IQR(variance))

