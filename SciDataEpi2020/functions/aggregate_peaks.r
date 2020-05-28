#### This is a poorly coded but perfectly functional (to the best of my knowledge) piece of code.

# This takes a GRanges object and searchs (locally) combining close bins with the previous by:
# keeping the smallest start bin position
# keeping the largest end bin position 
# either keeping the max score, or taking the average of all scores.
# defaults to 1000 DNA bases for search.

aggregate.peaks <- function(x,search.distance=1000,score.type="max"){
	peaks2keep <- data.frame(
				"chr"    = factor(seqnames(x)),
				"start"  = start(x),
				"end"    = end(x),
				"mcols"  =  mcols(x),
				"strand" = strand(x)
				)
	###################################################### MAX SCORING ######################################################
	if(score.type=="max"){
		for(i in 2:length(peaks2keep[,1])){
		# This checks to make sure that there is a peak to look at, then grabs the next peak ONLY if it's within search distance.
			while(if(i <= length(peaks2keep[,1])){abs(peaks2keep[i-1,"end"] - peaks2keep[i,"end"]) < search.distance}else{FALSE}){
			# find peak start of the two bins and keep the lowest
				peaks2keep[i-1,"start"] <- min(c(peaks2keep[i-1,"start"],peaks2keep[i,"start"]))
			# find peak end of the two bins and keep the highest
				peaks2keep[i-1,"end"]   <- max(c(peaks2keep[i-1,"end"],peaks2keep[i,"end"]))
			# find peak score of the two bins and keep the highest
			    peaks2keep[i-1,"mcols"] <- max(c(peaks2keep[i-1,"mcols"],peaks2keep[i,"mcols"]))    # for MAX scoring
			    
			# remove our bin that we've just extracted info from.
				peaks2keep <- peaks2keep[-i,]
			}
		}
	} else {
	###################################################### MEAN SCORING ######################################################
	for(i in 2:length(peaks2keep[,1])){
	# we need a counter for our average, resets every i.
		counter=1
		while(if(i <= length(peaks2keep[,1])){abs(peaks2keep[i-1,"end"] - peaks2keep[i,"end"]) < search.distance}else{FALSE}){
		# find peak start of the two bins and keep the lowest
			peaks2keep[i-1,"start"] <- min(c(peaks2keep[i-1,"start"],peaks2keep[i,"start"]))
		# find peak end of the two bins and keep the highest
			peaks2keep[i-1,"end"]   <- max(c(peaks2keep[i-1,"end"],peaks2keep[i,"end"]))
		# find peak score of the two bins and sum (we divide later after we calculate the counter
			peaks2keep[i-1,"mcols"] <- sum(peaks2keep[i-1,"mcols"],peaks2keep[i,"mcols"])         # for MEAN scoring
			# add one to counter which counts how many bins have been summed together
				counter = counter+1
		# remove our bin that we've just extracted info from.
				peaks2keep <- peaks2keep[-i,]
				head(peaks2keep,i)
			}
		# divide by the counter!
			if(i-1 <= length(peaks2keep[,1])) {peaks2keep[i-1,"mcols"] <- peaks2keep[i-1,"mcols"]/counter}
		}
	}
	###################################################### MAKE OUTPUT ######################################################
	output <-  GRanges(seqnames  = peaks2keep$chr, 
                        ranges   = IRanges(start=peaks2keep$start, end=peaks2keep$end),
                        strand   = peaks2keep$strand,
                        mcols    = peaks2keep$mcols
                        )
	return(output)
}


