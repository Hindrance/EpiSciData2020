	partition.samples <- function(classes, n){
 # find the minimal?
 
#	    n2sub <-  min(table(classes))
#	    sub.vector <- unlist(lapply(levels(classes2rf), function(k) {
#      sample(which(classes2rf==k), n2sub, replace=F)
#      }))
      
			output <- numeric(length(classes))
			train.vec <- unlist(lapply(1:length(unique(classes)), function(c) {sample(which(classes==unique(classes)[c]),(n*length(classes))/length(unique(classes)))}))
			output[train.vec] <- 1
			output[-train.vec] <- 2
			return(output)
			}
