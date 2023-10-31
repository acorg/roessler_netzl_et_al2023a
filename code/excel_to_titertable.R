# format titer table from excel sheet
rm(list = ls())
library(tidyverse)
library(stringr)
library(dplyr)

set_threshold <- function(tab, thresh) {
  tab[as.numeric(tab) < as.numeric(thresh)] <- paste0("<", thresh)
  tab[is.na(tab)] <- "*"
  
  return(tab)
}


# =============================== Omicron I sheet, published in NEJM
kimpel <- readxl::read_excel("./data/titer_data/230616_human vs hamster_update III.xlsx", range = "C2:O61", sheet = 1)
kimpel <- as.data.frame(kimpel)
kimpel[3,!is.na(kimpel[2,])] <- kimpel[2, !is.na(kimpel[2,])]
colnames(kimpel) <- kimpel[3,]
kimpel <- kimpel %>%
  filter(!is.na(`nAB ID`)) %>%
  filter(`nAB ID` != "nAB ID")

# get serum info
kimpel %>%
  mutate(inf_variant = str_extract(Sample, "[^ ]+"),
         inf_variant = gsub("first", "Anc. virus", inf_variant),
         inf_variant = gsub("Delta", "delta", inf_variant),
         sr_group = paste(inf_variant, "conv.", sep = " "), 
         sr_name = paste(sr_group, Species, `nAB ID`, Sample, sep = "_")) -> kimpel

# do human table
kimpel %>%
  filter(Species == "human") %>%
  select(D614G:sr_name) %>%
  select(!inf_variant:sr_group) %>%
  tibble::column_to_rownames("sr_name")-> human

human <- t(human)

human[as.numeric(human) < 16] <- "<16"
human[human == "n.a."] <- "*"

# do hamster
kimpel %>%
  filter(Species == "Hamster") %>%
  select(D614G:sr_name) %>%
  select(!inf_variant:sr_group) %>%
  tibble::column_to_rownames("sr_name")-> hamster

hamster <- t(hamster)
hamster[as.numeric(hamster) < 16] <- "<16"

# write tables
write.csv(human, "./data/titer_data/human_titer_table_<16.csv")
write.csv(hamster, "./data/titer_data/hamster_titer_table_<16.csv")

