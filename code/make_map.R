# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

# load data
mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 


rownames(mapColors) <- gsub("B.1.1.7", "Alpha", rownames(mapColors))
rownames(mapColors) <- gsub("P.1.1", "Gamma", rownames(mapColors))

alignment_map <- read.acmap("./data/maps/alignment_map.ace")
agNames(alignment_map) <- gsub(" omicron", "", agNames(alignment_map))
agNames(alignment_map) <- gsub("B.1.1.7", "Alpha", agNames(alignment_map))
agNames(alignment_map)[agNames(alignment_map) == "B.1.351"] <- "Beta"
agNames(alignment_map)[agNames(alignment_map) == "B.1.617.2"] <- "Delta"
agNames(alignment_map)[agNames(alignment_map) == "P.1.1"] <- "Gamma"
agNames(alignment_map)[agNames(alignment_map) == "XBB.1.5"] <- "XBB.1.5.1"
mapColors["XBB.1.5.1", "Color"] <- agFill(alignment_map)[agNames(alignment_map) == "XBB.1.5.1"]

make_map_from_table <- function(path_to_table, alignment_map){
  
  # load titer table
  titer_table <- read.titerTable(path_to_table)
 
  map <- make.acmap(titer_table, number_of_dimensions = 2, number_of_optimizations = 1000,
                      options =  list(ignore_disconnected = TRUE))
  
  dilutionStepsize(map) <- 0
  
  agFill(map) <- mapColors[agNames(map),]
  agGroups(map) <- agNames(map)
  
  sr_groups <- unlist(lapply(srNames(map), function(x) str_split(x, "_")[[1]][1]))
  
  srGroups(map) <- factor(sr_groups, levels = srGroup_colors$SerumGroup)
  srOutline(map) <- mapColors[as.character(srGroups(map)),]
  
  map <- apply_style(map)
  map<- realignMap(map, alignment_map)
  
  agSize(map)[agNames(map) == "Alpha+E484K"] <- unique(agSize(map))-3 
  
  return(map)
  
}



# ------------------------------------------ Make maps --------------------------------------
human_map <- make_map_from_table("data/titer_data/human_titer_table_<16.csv", alignment_map)
hamster_map <- make_map_from_table("data/titer_data/hamster_titer_table_<16.csv", human_map)

hamster_map <- reflectMap(hamster_map)
hamster_map <- rotateMap(hamster_map, 20)

save.acmap(human_map, "./data/maps/human_map_<16.ace")
save.acmap(hamster_map, "./data/maps/hamster_map_<16.ace")

