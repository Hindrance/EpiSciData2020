# Remove low expression function
	remove.low.expression.genes <- function(dataset, threshold.fraction = 0.5, bar=F){
		# threshold.fraction is the accepted number of 0 expression samples in a gene.
		# dataset is formatted as: dataset[genes,samples]
		dataset <- as.matrix(dataset)
		sample.no <- ncol(dataset)
		gene.no <- nrow(dataset)
		#pb <- txtProgressBar(0, gene.no, style=3)

		expression.fractions <- sapply(1:gene.no, function(i) {
						#setTxtProgressBar(pb,i)
						length(which(dataset[i,] == 0)) / sample.no
						}
					)
		cat("\n")
		expression.fractions.accepted <- expression.fractions < threshold.fraction
		output <- dataset[expression.fractions.accepted,]
		output
	}
	

