
# grab overlapping sections of genome where the GRanges object overlaps a reference data base containing features.
# example  - use the enseble "biomaRt" packages to download a data.frame of    | chromosome | start pos | end pos | gene  |

# where GR2 is the feature genomic ranges of interest with the mcols being the feature vector
	genomicFeatures <- function(GR1, GR2){
		n = length(GR1)
		x <- as.character(mcols(GR2[,1])[subjectHits(findOverlaps(GR1, GR2)),1])
		y <- queryHits(findOverlaps(GR1, GR2))
		z <- sapply(1:n, function(i) {paste(x[y==i],collapse=" | ")})
		return(z)
	}
	
