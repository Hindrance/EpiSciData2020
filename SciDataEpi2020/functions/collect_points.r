# This function calculates the density of overlapping points in 2-dimensions... It's relatively inefficient but is many times more obvious than
# other density functions, 2d histograms etc and results in a list of x + y coordinates, a density for the coordinates, and a colour scale using
# Rbrewer.

collect.points <- function(x, y, bin.dist=c(0,0), colours, method="round", weights, log=F, alpha=F){
  col2plot <- colorRampPalette(colors=colours, alpha=alpha)
    if(method=="round"){
      coord.density.a <- paste(round(x,bin.dist[1]),round(y,bin.dist[2]),sep="vs")
    } else {
    # create bin boundaries
      x.bins <- seq(round(min(x)-bin.dist[1], 10/bin.dist[1]), round(max(x)+bin.dist[1], 10/bin.dist[1]), by=bin.dist[1])
      y.bins <- seq(round(min(y)-bin.dist[2], 10/bin.dist[2]), round(max(y)+bin.dist[2], 10/bin.dist[2]), by=bin.dist[2])
    # calculate middle point of each bin
      x.bins.mid = x.bins + bin.dist[1]/2
      y.bins.mid = y.bins + bin.dist[2]/2
    # round the values in x to the x bin middles  
      x.binned <- character()
      for(bin in x.bins.mid){
        hits <- which(abs(bin-x) <= bin.dist[1]/2)
        x.binned[hits] <- as.character(bin)
      }
    # round the values in y to the y bin middles    
      y.binned <- character()
      for(bin in y.bins.mid){
        hits <- which(abs(bin-y) <= bin.dist[2]/2)
        y.binned[hits] <- as.character(bin)
      }
    # paste these together (a nice hack...)
      coord.density.a <- paste(x.binned,y.binned,sep="vs")
    }
  
  # table of counts for each x y combination (fast hack!)
	  coord.density <- table(coord.density.a)
	  # AMENDED 9th OCT 2019  - This is a weighting edit which allows us to aggrgate weight scores for each density bin... For example, the mean expression of a gene :).
	    if(!missing(weights)){
	     # if(sum(weights == 0) != length(weights)) {weights=normalise(weights)}
	     # if(sum(is.na(weights)) == length(weights)) {weights[1:length(weights)] = 0}
	      for(collected.bin in names(coord.density)){
	        # we calculate the mean of weights across all samples in bin region.
	        #coord.density[collected.bin] = coord.density[collected.bin] * (sum(weights[which(coord.density.a == collected.bin)]) / coord.density[collected.bin])
	        coord.density[collected.bin] = (sum(weights[which(coord.density.a == collected.bin)]) / coord.density[collected.bin])
        }
	    }
  # split by the separator (vs) 
	  test <- unlist(strsplit(names(coord.density), split="vs"))
	# format output :), calculating colours too...
	  data2plot <- data.frame(as.numeric(test[seq(1,length(test)-1, by =2)]), 
					  as.numeric(test[seq(2,length(test), by =2)]),
					  as.numeric(coord.density),
					  #as.character(col2plot(101)[round(normalise(log(as.numeric(coord.density)))*100)+1]),
					  if(log==T){as.character(col2plot(101)[as.numeric(cut(log(as.numeric(coord.density)), 101))])} else {as.character(col2plot(101)[as.numeric(cut(as.numeric(coord.density), 101))])},
					  stringsAsFactors=F
					  )
	  return(data2plot)
	}



collect.points.3d <- function(x, y, z, bin.dist=c(0,0,0), colours, method="round", weights, log=F){
  col2plot <- colorRampPalette(colors=colours)
    if(method=="round"){
      coord.density.a <- paste(round(x,bin.dist[1]),round(y,bin.dist[2]),round(z,bin.dist[3]),sep="vs")
    } else {
    # create bin boundaries
      x.bins <- seq(round(min(x)-bin.dist[1], 10/bin.dist[1]), round(max(x)+bin.dist[1], 10/bin.dist[1]), by=bin.dist[1])
      y.bins <- seq(round(min(y)-bin.dist[2], 10/bin.dist[2]), round(max(y)+bin.dist[2], 10/bin.dist[2]), by=bin.dist[2])
      z.bins <- seq(round(min(z)-bin.dist[3], 10/bin.dist[3]), round(max(z)+bin.dist[3], 10/bin.dist[3]), by=bin.dist[3])
    # calculate middle point of each bin
      x.bins.mid = x.bins + bin.dist[1]/2
      y.bins.mid = y.bins + bin.dist[2]/2
      z.bins.mid = z.bins + bin.dist[3]/2
    # round the values in x to the x bin middles  
      x.binned <- character()
      for(bin in x.bins.mid){
        hits <- which(abs(bin-x) <= bin.dist[1]/2)
        x.binned[hits] <- as.character(bin)
      }
    # round the values in y to the y bin middles    
      y.binned <- character()
      for(bin in y.bins.mid){
        hits <- which(abs(bin-y) <= bin.dist[2]/2)
        y.binned[hits] <- as.character(bin)
      }
    # round the values in z to the z bin middles    
      z.binned <- character()
      for(bin in z.bins.mid){
        hits <- which(abs(bin-z) <= bin.dist[3]/2)
        z.binned[hits] <- as.character(bin)
      }
    # paste these together (a nice hack...)
      coord.density.a <- paste(x.binned, y.binned, z.binned, sep="vs")
    }
  
  # table of counts for each x y combination (fast hack!)
	  coord.density <- table(coord.density.a)
	   # AMENDED 9th OCT 2019  - This is a weighting edit which allows us to aggrgate weight scores for each density bin... For example, the mean expression of a gene :).
	    if(!missing(weights)){
	      if(sum(weights == 0) == length(weights)) {weights = weights + 1}
	      #weights=normalise(weights)
	      if(sum(is.na(weights)) == length(weights)) {weights = weights + 1}
	      for(collected.bin in names(coord.density)){
	        # we calculate the mean of weights across all samples in bin region.
	        #coord.density[collected.bin] = coord.density[collected.bin] * (sum(weights[which(coord.density.a == collected.bin)]) / coord.density[collected.bin])
	        coord.density[collected.bin] = (sum(weights[which(coord.density.a == collected.bin)]) / coord.density[collected.bin])
        }
	    }
	    
  # split by the separator (vs) 
	  test <- unlist(strsplit(names(coord.density), split="vs"))
	# format output :), calculating colours too...
	  data2plot <- data.frame(
					  "x" = as.numeric(test[seq(1,length(test)-2, by =3)]), 
					  "y" = as.numeric(test[seq(2,length(test)-1, by =3)]),
					  "z" = as.numeric(test[seq(3,length(test),   by =3)]),
					  as.numeric(coord.density),
					  if(log==T){as.character(col2plot(101)[as.numeric(cut(log(as.numeric(coord.density)), 101))])} else {as.character(col2plot(101)[as.numeric(cut(as.numeric(coord.density), 101))])},
					  stringsAsFactors=F
					  )
	  return(data2plot)
	}
	





## Temporary work for colouring divergent sets - It basically uses two separate RGB or CMYK channels, one for each source.
## This is relatively simple for plotting data sets, however with collect points, the colours are calculated on a square / area basis.
#  w1 <- normalise(colMeans(d2norm.counts["BNC1",which(d2norm.dvds==day),drop=F]))
#  w2 <- normalise(colMeans(d2norm.counts["TCF21",which(d2norm.dvds==day),drop=F]))
#  
#  w1 <- colMeans(d2counts.DE[linA.markers,which(d2norm.dvds==day)])
#  w2 <- colMeans(d2counts.DE[linB.markers,which(d2norm.dvds==day)])
#  
#  w1 <- colMeans(d2counts.DE["BNC1",which(d2norm.dvds==day),drop=F])+min(w1[w1!=0])*0.000001
#  w2 <- colMeans(d2counts.DE["TCF21",which(d2norm.dvds==day),drop=F])+min(w2[w2!=0])*0.000001
#  w3 <- w1+w2
#  
#  
#  
#  w1 <- sort(runif(100, 0,1)); w1 = w1/max(w1)
#  w2 <- sort(runif(100, 0,500), decreasing=T); w2 = w2/max(w2)
#  w3 <- w1+w2
#  
#  w4 <- sapply(1:length(w1),function(i){w1[i]/max(c(w1[i],w2[i]))})
#  w5 <- sapply(1:length(w2),function(i){w2[i]/max(c(w1[i],w2[i]))})
#  w6 <- w3/max(w3)
#  col1 <- rgb(0.5+w4/2, 0.5+w5/2, 0.8, w6)
#  col1 <- rgb(w4, w5, 0, 1)
#  
#  plot(1:100,rep(1, 100),pch=15, col=col1,cex=1)

#  col1 <- rgb(w2,w2,w1)
#  col1 <- cmy(1-w1,1-w2,1)
#  plot(timepointPCAs[[i]]$x[,1], timepointPCAs[[i]]$x[,2], col=moo[,4], pch=NA, bg = "black", xlab="PC1", ylab="PC2", cex=1.3, asp=1)
#  rect(par("usr")[1], par("usr")[3],par("usr")[2],par("usr")[4],col=rgb(0,0,0))
#  #rect(par("usr")[1], par("usr")[3],par("usr")[2],par("usr")[4],col=cmy(0,0.1,0))
#  #points(timepointPCAs[[i]]$x[,1], timepointPCAs[[i]]$x[,2], bg=col1, pch=21, cex=1.3)
#  points(timepointPCAs[[i]]$x[,1], timepointPCAs[[i]]$x[,2], col=col1, pch=16, cex=1.3)
#  textplot()

#  w.mutual <- which(w1 > 0 & w2 > 0)
#  points(timepointPCAs[[i]]$x[w.mutual,1], timepointPCAs[[i]]$x[w.mutual,2], pch=1, col="white")
#  points(timepointPCAs[[i]]$x[w.mutual,1], timepointPCAs[[i]]$x[w.mutual,2], col=col1[w.mutual], pch=16, cex=1.3)
#  textplot(x=timepointPCAs[[i]]$x[w.mutual,1], y=timepointPCAs[[i]]$x[w.mutual,2], words=1:length(w[,1]), new=F, col="white")
#  
#  
#  w <- data.frame(w1[w1 > 0 & w2 > 0], w2[w1 > 0 & w2 > 0])
#  
#  points(timepointPCAs[[i]]$x[,1], timepointPCAs[[i]]$x[,2], col=gene.colour(d2norm.counts[,which(d2norm.dvds==day)], "TCF21"), pch=16)
#  points(timepointPCAs[[i]]$x[,1], timepointPCAs[[i]]$x[,2], col=gene.colour(d2norm.counts[,which(d2norm.dvds==day)], "BNC1"), pch=16)
#  
#  
#  points(c(1,2,3),c(1,1,1),pch=15, col=col1,cex=5)
#  
#  points(c(1,2,3),c(1,1,1),pch=15, col=rgb(0,0,1,0.5),cex=5)
#  points(c(1,2,3),c(1,1,1),pch=15, col=c(rgb(1,1,0,0.5), rgb(1,1,0,0.25), rgb(1,1,0,0)),cex=5)
#  
#  
#  w1 <- sort(runif(100, 0,1)); w1 = w1/max(w1)
#  w2 <- sort(runif(100, 0,500), decreasing=T); w2 = w2/max(w2)
#  w3 <- w1+w2
#  
#  w3 = (w1 - w2)/2+0.5
#  w4 = 0.5+normalise((w1 + w2)/2)/2
#  wcols <- colorRampPalette(colors=c("red", "yellow", "green"))(101) 

#  col1 = wcols[round(w3,2)*100+1]
#  
#  col1 <- rgb(w1,w2,0)
#  plot(1:100,rep(1, 100),pch=15, col=col1, cex=1)
#  points(1:100,rep(1, 100),pch=15, col=rgb(1,1,1,1-w4), cex=1)
#  









# OLD
## add one for 3d
#	collect.points.3d <- function(x,y,z,bin.dist=c(0,0,0),colours){
#	col2plot <- colorRampPalette(colors=colours)
#	coord.density <- paste(round(x,bin.dist[1]),round(y,bin.dist[2]),round(z,bin.dist[3]),sep="vs")
#	coord.density <- table(coord.density)
#	test <- unlist(strsplit(names(coord.density), split="vs"))
#	data2plot <- data.frame(
#					"x" = as.numeric(test[seq(1,length(test)-2, by =3)]), 
#					"y" = as.numeric(test[seq(2,length(test)-1, by =3)]),
#					"z" = as.numeric(test[seq(3,length(test),   by =3)]),
#					as.numeric(coord.density),
#					as.character(col2plot(101)[round(normalise(log(as.numeric(coord.density)))*100)+1]),
#					stringsAsFactors=F
#					)
#	return(data2plot)
#	}





















# OLD 
#	collect.points <- function(x,y,bin.dist=c(0,0),colours){
#	col2plot <- colorRampPalette(colors=colours)
#	coord.density <- paste(round(x,bin.dist[1]),round(y,bin.dist[2]),sep="vs")
#	coord.density <- table(coord.density)
#	test <- unlist(strsplit(names(coord.density), split="vs"))
#	data2plot <- data.frame(as.numeric(test[seq(1,length(test)-1, by =2)]), 
#					as.numeric(test[seq(2,length(test), by =2)]),
#					as.numeric(coord.density),
#					as.character(col2plot(101)[round(normalise(log(as.numeric(coord.density)))*100)+1]),
#					stringsAsFactors=F
#					)
#	return(data2plot)
#	}

