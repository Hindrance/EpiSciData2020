#geom.mean3 <- function(x, na.rm=TRUE){
#  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
#}

#geom.mean2 <- function(x) {
#	prod(x)^(1/length(x))
#}

## ignoring prod() flow issues that arise with geom.mean2
#geom.mean4 <- function(x) {
#              if (all(x == 0)) {
#                0
#              }
#              else {
#                exp(sum(log(x[x > 0]))/length(x))
#              }
#              
#            }
#            
#      

     
# THE reference function in DESeq
# https://rdrr.io/bioc/DESeq2/src/R/methods.R
geom.mean<- function(x,countMatrix=F){
	if(countMatrix == T) {
		return(exp(rowMeans(log(x))))
	} else {
		return(exp(mean(log(x))))
	}
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
