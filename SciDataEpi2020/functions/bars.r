
  bars <- function(data, labels=c(), error.data, label.cov = 0.66, bar.width, col, ylab="", xlab="", cex=1, ...){
    L <- length(data[,1])
    J <- length(data[1,])
      label.cov=0.5*label.cov
      
    bar.positions <- seq(0, label.cov, length=J) - mean(seq(0, label.cov, length=J))
    if(missing(bar.width)){bar.width = (label.cov / (J-1))}
    bp1 <- seq(-label.cov, label.cov-bar.width, length=J)
    bp2 <- seq(-label.cov+bar.width, label.cov, length=J)
  
    plot(seq(1-label.cov,L+label.cov,length=L),seq(0,max(data)*1.1,length=L), pch=NA, xaxt="n", ylab=ylab, xlab=xlab, cex.axis=cex, cex.lab=cex, ...)
    for(i in 1:L){
      for(j in 1:J){
        rect(i+bp1[j], 0, i+bp2[j], data[i,j], col=col[j],...)
      }
    }
    axis(1, at=(1:L), labels=rownames(data), cex.axis= cex, ...)
    
    if(!missing(error.data)) {
      for(i in 1:L){
        for(j in 1:J){
          arrows(mean(c(i+bp1[j],i+bp2[j])), data[i,j], mean(c(i+bp1[j],i+bp2[j])), data[i,j]+error.data[i,j], angle=90, length=0.1)
         #arrows(mean(c(i,i+bp1[j])), data[i,j], mean(c(i,i+bp1[j])), data[i,j]-error.data[i,j])
          #arrows(i+bp1[j], 0, i+bp2[j], data[i,j], col=col[j],...)
        }
      }
    }  
  }


