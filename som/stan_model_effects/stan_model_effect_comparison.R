rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggsci)
# Load data
sr_effect <- readRDS("data/titer_data/merged_map_stan_sampled_sr_imputed.rds") %>%
  select(source, sr_name, dataset_magnitude, sr_effect) %>%
  unique() %>%
  mutate(Model = "Serum reactivity")
org_effect <- readRDS("data/titer_data/merged_map_stan_sampled_magnitude_imputed.rds") %>%
  select(source, sr_name, dataset_magnitude, sr_effect) %>%
  unique() %>%
  mutate(Model = "Organism reactivity")
both <- readRDS("data/titer_data/merged_map_stan_sampled_both_imputed.rds") %>%
  select(source, sr_name, dataset_magnitude, sr_effect) %>%
  unique() %>%
  mutate(Model = "Serum & organism reactivity")


sr_org <- sr_effect %>%
  mutate(dataset_magnitude = org_effect$dataset_magnitude,
         sr_effect = sr_effect - dataset_magnitude,
         Model = "Serum with organism reactivity")

pal_npg(palette = c("nrc"), alpha = 1)
# combine all
combo <- rbind(sr_effect, org_effect, both) %>% #, sr_org) %>%
  mutate(sr_name = gsub("Anc. virus conv._|BA.1 conv._|BA.5 conv._|delta conv._|D614G conv._", "", sr_name),
         sr_name = gsub(" 1| 2| 3| 4|convalescent", " conv.", sr_name),
         sr_name = gsub("Hamster Sera \\(D21pi\\) Animal", "", sr_name),
         sr_name = gsub("Hamster Sera \\(D26pi\\) Animal", "", sr_name),
         sr_name = gsub("human", "Human", sr_name),
         sr_name = gsub("Delta", "delta", sr_name),
         sr_name = gsub("D614G|first wave", "Anc. virus", sr_name)
         )

combo$sr_name <- factor(combo$sr_name, levels = c(
  'Human_218_Anc. virus  conv.', 'Human_220_Anc. virus  conv.', 'Human_224_Anc. virus  conv.', 'Human_260_Anc. virus  conv.', 
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
  'Hamster_H2_BA.5  conv.', 'Hamster_H3_BA.5  conv.', 'Hamster_H4_BA.5  conv.'
))

# plot organism magnitude
combo %>%
  select(source, Model, dataset_magnitude) %>%
  unique() %>%
  ggplot(aes(x = source, y = dataset_magnitude, fill = Model)) + 
  geom_col(position = position_dodge(width = 0.7), color = "white") + 
  theme_bw() + 
  geom_hline(yintercept = 0, color = "grey20") +
  geom_point(position = position_dodge(width = 0.7), shape = 21, size = 3, color = "white")+
  scale_fill_npg() +
  scale_y_continuous(name = "Organism reactivity",
                     label = function(x) round(2^x, 2),
                     limits = c(-2.5,0.7),
                     breaks = c(-5:8)) + 
  xlab("Organism") -> p

ggsave("som/stan_model_effects/organism_effects.png", plot = p, dpi = 300, width = 5, height = 4)

# plot sr effects
combo %>%
  ggplot(aes(x = sr_name, y = sr_effect, color = Model, group = Model)) + 
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(name = "Serum reactivity",
                     label = function(x) round(2^x, 4),
                     limits = c(-12, 5),
                     breaks = c(-12:5)) + 
  xlab("Serum") + 
  scale_color_npg() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top") -> p

ggsave("som/stan_model_effects/sr_effects.png", plot = p, dpi = 300, width = 8, height = 4.7)

both %>%
  select(source, dataset_magnitude) %>%
  unique() %>%
  pivot_wider(names_from = source, values_from = dataset_magnitude) %>%
  View()
  mutate(org_reactivity_diff = 2^(Hamster-Human))

