
# Outlier expression detection??
	remove.outlier.mean.genes <- function(dataset, threshold) {
	  dataset <- as.matrix(dataset)
		dataset.log <- log(rowMeans(dataset))
		if(missing(threshold)){threshold = (mean(dataset.log)+2*sd(dataset.log))}
		cut.logX.pos <- mean(dataset.log) + threshold
		cut.logX.neg <- mean(dataset.log) - threshold
		x.2keep <- c(as.numeric(which(dataset.log < cut.logX.pos & dataset.log > cut.logX.neg)))
		output <- dataset[x.2keep,]
		output
	}		
	
	

