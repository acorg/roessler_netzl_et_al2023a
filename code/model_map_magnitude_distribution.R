# Setup workspace
rm(list = ls())
library(labbook)
library(patchwork)
library(Racmacs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(modeest)
source("functions/long_map_info.R")

set.seed(100)

# Function to model the map
model_map <- function(map, 
                      ag_mean_prior_mean = 3, # was 6
                      ag_mean_prior_sigma = 20,
                      sigma_prior_alpha = 2,
                      sigma_prior_beta = 5,
                      sigma_prior_sr_effect = 10,
                      sigma_prior_dataset_bias = 0.4, # was 0.4,
                      sigma_prior_dataset_magnitude = 2,
                      fitted_dataset_magnitude = NULL) {
  
  # Fetch data
  titer_matrix <- titerTable(map)
  ag_name_matrix <- Racmacs:::agNameMatrix(map)
  sr_extra_matrix <- srExtraMatrix(map)
  sr_name_matrix <- Racmacs:::srNameMatrix(map)
  sr_group_matrix <- Racmacs:::srGroupMatrix(map)
  
  # Vectorise data
  titers <- as.vector(titer_matrix)
  ag_name <- as.vector(ag_name_matrix)
  sr_name <- as.vector(sr_name_matrix)
  sr_extra <- as.vector(sr_extra_matrix)
  sr_group <- as.vector(sr_group_matrix)
  
  # Remove NA values
  na_titers <- titers == "*" | titers == "."
  titers <- titers[!na_titers]
  ag_name <- ag_name[!na_titers]
  sr_name <- sr_name[!na_titers]
  sr_extra <- sr_extra[!na_titers]
  sr_group <- sr_group[!na_titers]
  
  # Factor variables
  ag_name <- factor(ag_name)
  sr_name <- factor(sr_name)
  sr_extra <- factor(sr_extra)
  sr_group <- factor(sr_group)
  
  # Calculate titer limits
  titer_lims <- titertools:::calc_titer_lims(titers, dilution_stepsize = dilutionStepsize(map))
  
  # Build model
  model <- cmdstanr::cmdstan_model("models/titer_comparison_magnitude.stan")
  
  
  if(is.null(fitted_dataset_magnitude)){
    fitted_dataset_magnitude <- rep(0, nlevels(sr_extra))
  } 
  
  # Set input data
  stan_data <- list(
    N = length(titers),
    N_ags = nlevels(ag_name),
    N_srs = nlevels(sr_name),
    N_sr_groups = nlevels(sr_group),
    N_datasets = nlevels(sr_extra),
    upper_lims = titer_lims$max_titers,
    lower_lims = titer_lims$min_titers,
    ags = as.numeric(ag_name),
    srs = as.numeric(sr_name),
    sr_groups = as.numeric(sr_group),
    datasets = as.numeric(sr_extra),
    ag_mean_prior_mean = ag_mean_prior_mean,
    ag_mean_prior_sigma = ag_mean_prior_sigma,
    sigma_prior_alpha = sigma_prior_alpha,
    sigma_prior_beta = sigma_prior_beta,
    sigma_prior_sr_effect = sigma_prior_sr_effect,
    sigma_prior_dataset_bias =sigma_prior_dataset_bias,
    sigma_prior_dataset_magnitude = sigma_prior_dataset_magnitude, 
    dataset_magnitude_mean = fitted_dataset_magnitude  
  )
  
  if(is.null(fitted_dataset_magnitude)){
    fitted_dataset_magnitude <- rep(0, nlevels(sr_extra))
  } 
 
  ag_means <- matrix(NA, nlevels(ag_name), nlevels(sr_group))
  ag_level_matrix <- matrix(levels(ag_name), nlevels(ag_name), nlevels(sr_group), byrow = F)
  sr_group_level_matrix <- matrix(levels(sr_group), nlevels(ag_name), nlevels(sr_group), byrow = T)
  ag_means[] <- vapply(seq_along(ag_means), \(i){
    mean(titer_lims$log_titers[ag_name == ag_level_matrix[i] & sr_group == sr_group_level_matrix[i]])
  }, numeric(1))
  ag_means[is.na(ag_means)] <- stan_data$ag_mean_prior_mean
  
  stan_init <- list(
    ag_means = ag_means,
    sr_effects = rep(0, nlevels(sr_name)),
    dataset_bias = array(0, c(nlevels(ag_name), nlevels(sr_group), nlevels(sr_extra))),
    dataset_magnitude = fitted_dataset_magnitude,
    sigma_error = matrix(1, nlevels(sr_extra), nlevels(sr_group))
  )
  
  # Optimize model
  model_optim <- model$optimize(
    data = stan_data,
    seed = 100,
    init = list(
      stan_init
    )
  )
  
  
  fit <- model$sample(
    data = stan_data,
    seed = 100,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    chains = 4,
    parallel_chains = 4,
    refresh = 500 # print update every 500 iters
  )
  
 # return(list("fit" = fit, "optim" = model_optim, "model" = model))
  model_optimized <- fit
    # Get generated quantities
    model_pars <- model_optimized$summary()[-1,]
    model_par_list <- as.list(model_pars$mean)
    names(model_par_list) <- model_pars$variable
    
    imputed_logtiters <- model$generate_quantities(
      fitted_params = do.call(posterior::draws_matrix, model_par_list),
      data = stan_data
    )$summary()
    
    imputed_logtiters |> 
      mutate(
        sr_name = sr_name,
        ag_name = ag_name
      ) |> 
      rename(
        imputed_logtiter = mean
      ) |> 
      select(
        sr_name, 
        ag_name,
        imputed_logtiter
      ) -> imputed_logtiters
    
    # Organize results
    sr_effects <- model_optimized$summary("sr_effects") |> mutate(sr_name = levels(sr_name)) |> select(sr_name, mean) |> rename(sr_effect = mean)
    dataset_magnitude <- model_optimized$summary("dataset_magnitude") |> mutate(source = levels(sr_extra)) |> select(source, mean) |> rename(dataset_magnitude = mean)
    
    ag_gmts <- model_optimized$summary("ag_means") |>
      mutate(
        indices = gsub("^.*\\[", "[", variable),
        indices = gsub("\\[", "", indices),
        indices = gsub("\\]", "", indices)
      ) |> 
      separate(
        indices,
        sep = ",",
        into = c("ag", "group")
      ) |> 
      mutate(
        ag_name = as.character(levels(ag_name)[as.numeric(ag)]),
        sr_group = as.character(levels(sr_group)[as.numeric(group)])
      ) |> 
      select(
        ag_name, sr_group, mean
      ) |> rename(
        ag_mean = mean
      )
    
    dataset_error <- model_optimized$summary("sigma_error") |> 
      mutate(
        indices = gsub("^.*\\[", "[", variable),
        indices = gsub("\\[", "", indices),
        indices = gsub("\\]", "", indices)
      ) |> 
      separate(
        indices,
        sep = ",",
        into = c("dataset", "group")
      ) |> 
      mutate(
        source = as.character(levels(sr_extra)[as.numeric(dataset)]),
        sr_group = as.character(levels(sr_group)[as.numeric(group)])
      ) |> 
      select(
        source, sr_group, mean
      ) |> 
      rename(
        dataset_error = mean
      )
    
    dataset_bias <- model_optimized$summary("dataset_bias") |> 
      mutate(
        indices = gsub("^.*\\[", "[", variable),
        indices = gsub("\\[", "", indices),
        indices = gsub("\\]", "", indices)
      ) |> 
      separate(
        indices,
        sep = ",",
        into = c("ag", "group", "dataset")
      ) |> 
      mutate(
        ag_name = as.character(levels(ag_name)[as.numeric(ag)]),
        sr_group = as.character(levels(sr_group)[as.numeric(group)]),
        source = as.character(levels(sr_extra)[as.numeric(dataset)])
      ) |> 
      select(
        source, sr_group, ag_name, mean
      ) |> 
      rename(
        dataset_sr_group_bias = mean
      )
    
    # Get long version of the map data
    mapdata <- long_map_info(map) |> rename(source = sr_extra)
    
    # Bind in estimates
    mapdata |> 
      left_join(sr_effects, by = "sr_name") |> 
      left_join(dataset_magnitude, by = "source") |> 
      left_join(dataset_bias, by = c("ag_name", "source", "sr_group")) |> 
      left_join(imputed_logtiters, by = c("ag_name", "sr_name")) |> 
      left_join(ag_gmts, by = c("ag_name", "sr_group")) -> mapdata
    
    return(list("mapdata" = mapdata, "fit" = fit, "optim" = model_optim))
 
  # 
  
  
}


hamster_map <- read.acmap("./data/maps/hamster_map_<16.ace")
human_map <- read.acmap("./data/maps/human_map_<16.ace")

srExtra(hamster_map) <- rep("Hamster", length(srExtra(hamster_map)))
srExtra(human_map) <- rep("Human", length(srExtra(human_map)))
merged_table <- mergeMaps(list(human_map, hamster_map))
merged_table_adjusted <- cbind(titerTable(human_map), titerTable(hamster_map))
merged_map <- make.acmap(merged_table_adjusted, number_of_dimensions = 2, number_of_optimizations = 0, options = list(ignore_disconnected = TRUE))
srExtra(merged_map) <- srExtra(merged_table)
srGroups(merged_map) <- srGroups(merged_table)
srGroups(merged_map) <- gsub("D614G", "Anc. virus", srGroups(merged_map))

dilutionStepsize(merged_map) <- 0

# Model first magnitude effects only
magnitude_og <- model_map(merged_map, sigma_prior_sr_effect = 1e-3, sigma_prior_dataset_bias = 1e-3, sigma_prior_dataset_magnitude = 10) 

saveRDS(
  object = magnitude_og$mapdata,
  "data/titer_data/merged_map_stan_sampled_magnitude_imputed.rds"
)

sr_og <- model_map(merged_map, sigma_prior_sr_effect = 10, sigma_prior_dataset_bias = 1e-3, sigma_prior_dataset_magnitude = 1e-3) 

saveRDS(
  object = sr_og$mapdata,
  "data/titer_data/merged_map_stan_sampled_sr_imputed.rds"
)


magnitude_both <- model_map(merged_map, sigma_prior_sr_effect = 10, sigma_prior_dataset_bias = 1e-3, sigma_prior_dataset_magnitude = 10) 

saveRDS(
  object = magnitude_both$mapdata,
  "data/titer_data/merged_map_stan_sampled_both_imputed.rds"
)

#------------------ Do data formatting for plot
titer_matrix <- titerTable(merged_map)
titers <- as.vector(titer_matrix)

# get reference infos
sr_extra_matrix <- srExtraMatrix(merged_map)
sr_extra <- as.vector(sr_extra_matrix)

# ag infos
ag_name_matrix <- Racmacs:::agNameMatrix(merged_map)
sr_name_matrix <- Racmacs:::srNameMatrix(merged_map)
sr_group_matrix <- Racmacs:::srGroupMatrix(merged_map)

ag_name <- as.vector(ag_name_matrix)
sr_name <- as.vector(sr_name_matrix)
sr_group <- as.vector(sr_group_matrix)

# Remove NA values
na_titers <- titers == "*" | titers == "."
titers <- titers[!na_titers]
ag_name <- ag_name[!na_titers]
sr_name <- sr_name[!na_titers]
sr_extra <- sr_extra[!na_titers]
sr_group <- sr_group[!na_titers]


# Factor variables
ag_name <- factor(ag_name)
sr_name <- factor(sr_name)
sr_group <- factor(sr_group)
sr_extra <- factor(sr_extra)

sr_extra_levels <- levels(sr_extra)
sr_group_levels <- levels(sr_group)
ag_means <- matrix(NA, nlevels(ag_name), nlevels(sr_group))

ag_level_matrix <- matrix(levels(ag_name), nlevels(ag_name), nlevels(sr_group), byrow = F)
sr_group_level_matrix <- matrix(levels(sr_group), nlevels(ag_name), nlevels(sr_group), byrow = T)

for(r in 1:nrow(ag_means)){
  for(c in 1:ncol(ag_means)){
      ag_means[r,c] <- paste0(ag_level_matrix[r, c], ", ", sr_group_level_matrix[r, c])
  }
}

sr_effects <- levels(sr_name)
dataset_magnitude <- sr_extra_levels


plot_distribution_and_optim <- function(models, target_var = "dataset_magnitude|sr_effect|ag_means"){
  
  model_optimized <- models$fit
  mean_val <- model_optimized$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  ) %>%
    mutate(keep = grepl(target_var, variable)) %>%
    filter(keep) %>%
    pivot_longer(cols = c("mean", "median"), names_to = "Data", values_to = "estimate") %>%
    select(keep, Data, estimate, variable)
  
  optim_val <- models$optim$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
  ) %>%
    mutate(keep = grepl(target_var, variable)) %>%
    filter(keep) %>%
    mutate("Data" = "optimization") %>%
    rbind(., mean_val)
  
  draws_df <- model_optimized$draws(format = "df")
 
  sub_df <- draws_df[, grepl(target_var, colnames(draws_df))]
  sub_df %>%
    pivot_longer(cols = colnames(sub_df), names_to = "variable", values_to = "estimate") -> sub_df_long
  
  sub_df_long$variable <- sapply(sub_df_long$variable, function(x) eval(parse(text = x)))
 
  optim_val$variable <- sapply(optim_val$variable, function(x) eval(parse(text = x)))
  
  x_var_labels <- c("sr_effect" = "Serum",
                    "dataset_magnitude" = "Organism",
                    "ag_means" = "Variant, serum group")
  
  sr_levels <- c('Human_218_Anc. virus  conv.', 'Human_220_Anc. virus  conv.', 'Human_224_Anc. virus  conv.', 'Human_260_Anc. virus  conv.', 
                 'Human_262_Anc. virus  conv.', 'Human_278_Anc. virus  conv.', 'Human_289_Anc. virus  conv.', 'Human_280_Anc. virus  conv.', 
                 'Human_298_Anc. virus  conv.', 'Human_299_Anc. virus  conv.', 
                 'Human_G21_delta  conv.', 'Human_G22_delta  conv.', 'Human_G23_delta  conv.', 'Human_G24_delta  conv.', 'Human_G25_delta  conv.', 
                 'Human_G26_delta  conv.', 'Human_G27_delta  conv.', 
                 'Human_G291_BA.1   conv.', 'Human_G347_BA.1   conv.', 'Human_G348_BA.1   conv.', 'Human_G353_BA.1   conv.', 'Human_G354_BA.1   conv.', 
                 'Human_G355_BA.1   conv.', 'Human_G356_BA.1   conv.', 'Human_G359_BA.1   conv.', 'Human_G360_BA.1   conv.', 'Human_G362_BA.1   conv.', 
                 'Human_G363_BA.1   conv.', 'Human_G366_BA.1   conv.', 'Human_G369_BA.1   conv.', 'Human_G371_BA.1   conv.', 'Human_G381_BA.1   conv.', 
                 'Human_G393_BA.1   conv.', 'Human_G647_BA.1   conv.', 'Human_G650_BA.1   conv.', 
                 'Human_G790_BA.5  conv.', 'Human_G938_BA.5  conv.', 'Human_G939_BA.5  conv.', 
                 'Hamster_G990_Anc. virus  conv.', 'Hamster_G991_Anc. virus  conv.', 'Hamster_G992_Anc. virus  conv.', 'Hamster_G993_Anc. virus  conv.', 
                 'Hamster_G994_delta  conv.', 'Hamster_G995_delta  conv.', 'Hamster_G996_delta  conv.', 'Hamster_G997_delta  conv.', 
                 'Hamster_G998_BA.1  conv.', 'Hamster_G999_BA.1  conv.', 'Hamster_G1000_BA.1  conv.', 'Hamster_H1_BA.1  conv.', 
                 'Hamster_H2_BA.5  conv.', 'Hamster_H3_BA.5  conv.', 'Hamster_H4_BA.5  conv.')
  
  sub_df_long %>%
    mutate(variable = gsub("Anc. virus conv._|BA.1 conv._|BA.5 conv._|delta conv._|D614G conv._", "", variable),
           variable = gsub(" 1| 2| 3| 4|convalescent", " conv.", variable),
           variable = gsub("Hamster Sera \\(D21pi\\) Animal", "", variable),
           variable = gsub("Hamster Sera \\(D26pi\\) Animal", "", variable),
           variable = gsub("human", "Human", variable),
           variable = gsub("Delta", "delta", variable),
           variable = gsub("D614G|first wave", "Anc. virus", variable)
    ) -> sub_df_long
  
  optim_val %>%
    mutate(variable = gsub("Anc. virus conv._|BA.1 conv._|BA.5 conv._|delta conv._|D614G conv._", "", variable),
           variable = gsub(" 1| 2| 3| 4|convalescent", " conv.", variable),
           variable = gsub("Hamster Sera \\(D21pi\\) Animal", "", variable),
           variable = gsub("Hamster Sera \\(D26pi\\) Animal", "", variable),
           variable = gsub("human", "Human", variable),
           variable = gsub("Delta", "delta", variable),
           variable = gsub("D614G|first wave", "Anc. virus", variable)
    ) -> optim_val
  
  if(target_var == "sr_effect"){
    optim_val$variable <- factor(optim_val$variable, levels = sr_levels)
    sub_df_long$variable <- factor(sub_df_long$variable, levels = sr_levels)
  }
  
  
  sub_df_long %>%
    ggplot(aes(y = estimate, x = variable)) + 
    geom_violin() + 
    geom_point(data = optim_val, aes(color = Data)) + 
    theme_bw() +
    ylab("Estimate") + 
    xlab(x_var_labels[target_var]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))-> pl
  
  if(target_var == "sr_effect"){
    pl <- pl + 
    scale_y_continuous(name = "Serum reactivity",
                       label = function(x) round(2^x, 3),
                       limits = c(-10, 10))  
  } else if(target_var == "ag_means") {
    pl <- pl + 
      scale_y_continuous(name = "Geometric mean titer",
                         label = function(x) round(2^x*10),
                         limits = c(-10, 25))  
  } else {
    pl <- pl + 
      scale_y_continuous(name = "Organism reactivity",
                         label = function(x) round(2^x, 3),
                         limits = c(-15, 8))  
  }
  
  return(pl)
}


plot_all_distributions <- function(model){
  
  plot_list <- list()
  for(var in c("ag_means", "sr_effect", "dataset_magnitude")){
    plot_list[[var]] <- plot_distribution_and_optim(model, var)
  }
 
  wrapped <- wrap_plots(plot_list, nrow = 3) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')
  return(wrapped)
  
}

plots <- plot_all_distributions(magnitude_both)
ggsave("som/stan_model_effects/effect_distribution_both.png", plots, width = 12, height = 14, dpi = 300)

plots <- plot_all_distributions(magnitude_og)
ggsave("som/stan_model_effects/effect_distribution_magnitude.png", plots, width = 12, height = 14, dpi = 300)


plots <- plot_all_distributions(sr_og)
ggsave("som/stan_model_effects/effect_distribution_serum.png", plots, width = 12, height = 14, dpi = 300)

