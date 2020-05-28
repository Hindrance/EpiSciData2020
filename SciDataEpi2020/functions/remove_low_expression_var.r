
	
# remove low variance genes	
remove.low.variance.genes <- function(dataset, threshold.variance = 2){
		# threshold.fraction is the accepted number of 0 expression samples in a gene.
		# dataset is formatted as: dataset[genes,samples]
		dataset <- as.matrix(dataset)
		sample.no <- ncol(dataset)
		gene.no <- nrow(dataset)
		#pb <- txtProgressBar(0, gene.no, style=3)

		
		ev <- sapply(1:gene.no, function(i) {
			#setTxtProgressBar(pb,i)					
                        max(dataset[i,]) / min(dataset[i,])
			}
		)		
		cat("\n")
		variances.accepted <- which(ev >= threshold.variance)
		output <- dataset[variances.accepted,]
		output
	}


