rm(list = ls())
library(Racmacs)
set.seed(100)

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x


for(org in c("human", "hamster")){
  cat(org, " start")
  neut <- read.acmap(paste0('data/maps/', org, '_map_<16.ace'))
  
  
  # Titer and antigen noise
  neutBootTA <- bootstrapMap(
    neut,
    "bayesian",
    bootstrap_repeats = 500, # was 1000
    bootstrap_ags = TRUE,
    bootstrap_sr = TRUE,
    reoptimize = TRUE,
    optimizations_per_repeat = 1000,
    ag_noise_sd = 0.7,
    titer_noise_sd = 0.7,
    options = list(ignore_disconnected = TRUE)
  )
  
  save.acmap(neutBootTA, paste0("./som/bootstrapping/", org, "_neutBootTA_bayesian_unadj.ace"))
  
  cat("TA done")
  
}
  # do the blobs
  neutBootTA <- read.acmap(paste0("./som/bootstrapping/human_neutBootTA_bayesian_unadj.ace"))
  neutBootTABlobs <- bootstrapBlobs(neutBootTA, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
  neutBootT<- read.acmap(paste0("./som/bootstrapping/hamster_neutBootTA_bayesian_unadj.ace"))
  neutBootTBlobs <- bootstrapBlobs(neutBootT, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
  
  # Plot the figure
  png(paste0("./som/bootstrapping/both_unadj_bayesian-bootstrap.png"), width = 6, height = 3, units = 'in', res=300, pointsize = 18)
  layout(matrix(c(1, 2), ncol=2))
  par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
  plot(neutBootTABlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
  text(xlim_no_zoom[1]+0.6, ylim_no_zoom[2]-0.6, 'A', cex = 1.4)
  plot(neutBootTBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
  text(xlim_no_zoom[1]+0.6, ylim_no_zoom[2]-0.6, 'B', cex = 1.4)
  dev.off()
  

