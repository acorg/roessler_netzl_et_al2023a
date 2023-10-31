# labelled map 
# Setup plotting function
doplot <- function(map, xlims, ylims, show_labels = TRUE, ag_arrow_coords = NULL,
                   sr_arrow_coords = NULL) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  srSize(map) <- srSize(map) - 2
  agSize(map) <- agSize(map) - 2
  agSize(map)[agNames(map) == "Alpha+E484K"] <- agSize(map)[agNames(map) == "Alpha+E484K"]-3
  # Plot the regular map
  srOutlineWidth(map) <- 0.8
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.8)
  plot(map, xlim = xlims, 
       ylim =ylims, fill.alpha = 0.9,
       plot_stress = TRUE) ->p
  
  if(!is.null(ag_arrow_coords)){
    for(ar in ag_arrow_coords){
      arrows(x0 = ar[["x0"]], y0 = ar[["y0"]], x1 = ar[["x1"]], y1= ar[["y1"]], length = 0.1 )
    }
  }
  
  return(p)

}

