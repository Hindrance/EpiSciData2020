#   This function produces a histogram with a rectangle on either side of a 
#   bounds vector suggesting the regions that are being filtered out or selected against.


  hist.boundaries <- function(x, x.bounds, breaks = 50, colour = "red", new=T, ...){
      H <- hist(x, breaks = 50, plot=F)
      if(new==T) {plot(H, ylim=c(0,max(H$counts)*1.1), ...)}
      x.range = diff(range(x))
      plot.lim = x.range * 0.1
      rect(min(x) - plot.lim, 0, x.bounds[1],  max(H$counts), density=20, col=colour)
      rect(x.bounds[2], 0, max(x)+plot.lim, max(H$counts), density=20, col=colour)
      text(mean(x[x < x.bounds[1]]), max(H$counts)*1.02, sum(x < x.bounds[1]), xpd=NA, col=colour)
      text(mean(x[x > x.bounds[1] & x < x.bounds[2]]), max(H$counts)*1.08, sum(x > x.bounds[1] & x < x.bounds[2]), xpd=NA)
      text(mean(x[x > x.bounds[2]]), max(H$counts)*1.02, sum(x > x.bounds[2]), xpd=NA, col=colour) 
  }
