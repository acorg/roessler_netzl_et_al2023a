# Setup workspace
rm(list = ls())
library(labbook)
library(tidyverse)
library(patchwork)
library(Racmacs)
library(ggh4x)

source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")


# Function to plot results split by source
plot_sr_group_results_sub <- function(mapdata, ag_order = NULL, facet_by = "facet_col", ag_means_all = FALSE) {
  
  # Apply adjustments
  mapdata$source <- factor(mapdata$source)
  
  # Get gmt data
  ag_gmts <- mapdata |> distinct(ag_name, sr_group, ag_mean) 
  
  mapdata -> mapdata_subset
  
  ag_gmts |> 
    filter(
      ag_name %in% filter(mapdata_subset, !is.na(logtiter))$ag_name
    ) %>%
    mutate(facet_col = "Human") -> ag_gmts_sub
  
  ag_gmts_subset <- rbind(ag_gmts_sub, ag_gmts_sub %>%
                     mutate(facet_col = "Hamster"))
  
  if(ag_means_all){
    ag_gmts_subset <- rbind(ag_gmts_subset, ag_gmts_sub %>%
                              mutate(facet_col = "GMT"))
  }
  
  mapdata_subset |> 
    filter(
      ag_name %in% filter(mapdata_subset, !is.na(logtiter))$ag_name
    ) -> mapdata_subset
  
  # Set antigen order
  if(is.null(ag_order)){
    ag_order <- ag_gmts_subset |> arrange(-ag_mean) |> pluck("ag_name") 
  }
  
  # Do the plot
  mapdata_subset |>
    mutate(
      plot_type = "Raw titers"
    ) -> plotdata_raw
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - sr_effect ,
      plot_type = "Adj. for serum reactivity"
    ) -> plotdata_sr_effect
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - sr_effect - dataset_magnitude,
      plot_type = "Adj. for serum & organism reactivity"
    ) -> plotdata_sr_effect_magnitude
  
  mapdata_subset |>
    mutate(
      logtiter = imputed_logtiter - dataset_magnitude,
      plot_type =  "Adj. for organism reactivity"
    ) -> plotdata_magnitude
  
  plotdata <- bind_rows(
    plotdata_raw,
    plotdata_sr_effect,
    plotdata_sr_effect_magnitude,
    plotdata_magnitude
  ) |> 
    mutate(
      facet_col = source,
      plot_type = factor(
        plot_type,
        c(
          "Raw titers",
          "Adj. for serum reactivity",
          "Adj. for organism reactivity",
          "Adj. for serum & organism reactivity"
        )
      )
    )
  
  # calculate gmt now
  plotdata %>%
    group_by(source, plot_type, ag_name, sr_group) %>%
    mutate(logtiter = ifelse(titer == "<16", log2(8/10), logtiter),
            lower = Rmisc::CI(na.omit(logtiter))["lower"],
           upper = Rmisc::CI(na.omit(logtiter))["upper"],
           logtiter = mean(logtiter, na.rm =TRUE)) %>%
    ungroup() %>%
    mutate(facet_col = "GMT",
           sr_name = paste0("GMT", source)) %>%
    select(source, plot_type, ag_name, sr_group, sr_name, lower, upper, logtiter, facet_col) %>%
    unique() %>%
    mutate(alpha_val = "GMT")-> plotdata_source
  
  if(facet_by == "source"){
    plotdata_source <- rbind(plotdata_source, plotdata_source %>%
                               mutate(facet_col = source)) 
  }

  plotdata <- plyr::rbind.fill(plotdata %>%
                                 mutate(alpha_val = "Organism"), plotdata_source)
 
 plotdata %>%
  # mutate(alpha_val = ifelse(facet_col == "GMT", "GMT", "Organsim")) -> plotdata
  mutate(alpha_val = ifelse(grepl("GMT", sr_name), "GMT", "Organsim"),
         facet_col = factor(facet_col, levels = c("Human", "Hamster", "GMT")),
         sr_group = factor(sr_group, levels = c("Anc. virus conv.", "delta conv.", "BA.1 conv.", "BA.5 conv."))) -> plotdata 
 
 ag_gmts_subset <- 
   ag_gmts_subset %>%
   mutate(facet_col = factor(facet_col, levels = c("Human", "Hamster", "GMT")),
          sr_group = factor(sr_group, levels = c("Anc. virus conv.", "delta conv.", "BA.1 conv.", "BA.5 conv.")))
 
 plotdata |> 
    ggplot(
      aes(
        x = ag_name,
        y = logtiter,
        group = sr_name,
        color = source,
        alpha = alpha_val
      )
    ) + 
    geom_ribbon(
      aes(ymin = lower,
          ymax = upper,
          fill = source),
      color = NA,
      alpha = 0.3
    ) + 
    geom_line(
     # alpha = 0.4
    ) + 
    geom_line(
      data = plotdata |> filter(!is.na(logtiter)),
    #  alpha = 0.4,
      linetype = "dotted"
    ) + 
    scale_alpha_discrete(range = c(0.4, 1), 
                         limits = c("Organsim", "GMT")) + 
    geom_line(
      data = ag_gmts_subset,
      mapping = aes(
        y = ag_mean,
        group = "average"
      ),
      color = "black",
      alpha =1
    ) +
    annotate(
      "tile",
      x = unique(mapdata$ag_name),
      y = -2,
      height = 2*(2 + log2(16/10)),
      fill = "black",
      color = NA,
      alpha = 0.1
    ) +
    scale_x_discrete(
      limits = ag_order,
      expand = c(0, 0)
    ) +
    scale_y_continuous(labels = function(x) 2^x*10,
                       breaks = seq(-1,12,2)) +
    scale_color_manual(values = c("Hamster" = "#F8766D",
                                  "Human" = "#619CFF")) +
    scale_fill_manual(values = c("Hamster" = "#F8766D",
                                  "Human" = "#619CFF")) + 
    facet_nested(plot_type ~ sr_group + facet_col,
      drop = FALSE,
      labeller = label_wrap_gen(multi_line = TRUE),
      nest_line = element_line(color = "grey20")
    ) +
    coord_cartesian(
      ylim = c(-2, 12)
    ) +
    labs(
      x = "",
      y = "Titer"
    ) + 
    titerplot_theme() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.text = element_text(size = 9),
          strip.text.y= element_text(size = 9),
          strip.text.x= element_text(size = 9))-> gp1
  
  gp <- gp1
  
  return(gp)
}

figure_dir <- file.path("figures", "titerplots_stan")
if(!exists(figure_dir)){
  dir.create(figure_dir)
}


human_map <- read.acmap("./data/maps/human_map_<16.ace")


for(data in c("_sampled_both_imputed", "_sampled_sr_imputed","_sampled_magnitude_imputed")){ 
 
  mapdata <- readRDS(paste0("data/titer_data/merged_map_stan", data, ".rds"))

  wrapped_sub <- plot_sr_group_results_sub(mapdata, ag_order = agNames(human_map), facet_by = "facet_col", ag_means_all = TRUE)
  
  ggsave(file.path(figure_dir, paste0("facet_ag_order_titerplots_stan_gmt",data,"_", "facet_col", "_agGMT_", TRUE,".png")), plot = wrapped_sub, width = 18, height = 8)
  
 
}

