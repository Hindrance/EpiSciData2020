 # This function adds a gradient legend to the side of a plot. So far it only includes
 # the option of side 4 (right hand side).


 gradient.legend <- function(data, side=4, colours, scale=1, distance=0, label="", xpd=NA, cex=1, col.axis="black"){
 # some adjustments to ensure that it does not override current par() settings...
# cex.old = par()$cex
 #par(cex=cex)
 
  x.lim <- par("usr")[2]*((distance*0.08)+1.05)
  
    colour.scale <- rev(rgb(colorRamp(colours,alpha=T)(normalise(1:200))/255,alpha=(colorRamp(colours,alpha=T)(normalise(1:200))/255)[,4]))
    legend.x <- rep(x.lim,200)#*((distance*0.08)+1.05)
    legend.y <- seq(par("usr")[4]-(1-scale)*diff(c(par("usr")[3],par("usr")[4])), par("usr")[3]+(1-scale)*diff(c(par("usr")[3],par("usr")[4])), length=200)
#    axis(side, pos=x.lim, at=seq(par("usr")[4]-(1-scale)*diff(c(par("usr")[3],par("usr")[4])), par("usr")[3]+(1-scale)*diff(c(par("usr")[3],par("usr")[4])),length=10), labels=signif(seq(min(data), max(data), length=10),2),line=1.1+distance)
    axis(side, pos = x.lim*1.02, at=mean(par("usr")[3:4]), labels=label, tcl = 0, col.ticks="white") 
    axis(side, pos=x.lim, at=seq(par("usr")[4]-(1-scale)*diff(c(par("usr")[3],par("usr")[4])), par("usr")[3]+(1-scale)*diff(c(par("usr")[3],par("usr")[4])),length=10), labels=signif(seq(max(data), min(data), length=10),2), padj=-1, cex.axis=cex, col.axis=col.axis, col.ticks=col.axis, col=col.axis, col.lab=col.axis)
    points(legend.x, legend.y, pch=15, col=colour.scale, xpd=NA)
    
  # reset
# par()$cex = cex.old
  }
  

