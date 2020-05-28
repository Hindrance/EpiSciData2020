DESeq.normalise <- function(counts.matrix){

# Turn it into a matrix and calculate dimensions
	counts.matrix <- as.matrix(counts.matrix)
	n = length(counts.matrix[,1]) # rows are genes
	J = length(counts.matrix[1,]) # columns are samples

# calculate geometric mean of each gene across samples
	a <<- sapply(1:n, function(i) {geom.mean(counts.matrix[i,])}) 

# divide each gene's read counts in each sample by the population geometric mean where counts are above 0.
	scale.factor <<- sapply(1:J, function(j) {
		median((counts.matrix[,j]/a)[a>0])
		}
	)
	
# divide each sample by its scale factor
	b <- t(t(as.matrix(counts.matrix))/scale.factor)
	
	return(b)
	}
	
