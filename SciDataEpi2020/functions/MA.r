

MA <- function(data, size, outsize=T){
	if(missing(size)){size=3}
	data <- as.numeric(data)
	L <- length(data)
	output <- numeric(L)
	if(outsize==T){
		for(i in 1:L) {
			data.i <- (i-floor(size/2)):((i-1)+round(size/2))
				if(length(data.i) != size) {
					size.i = (length(which(data.i > 0 & data.i < L)) - length(which(!(data.i > 0 & data.i < L))))
					data.i <- (i-floor(size.i/2)):((i-1)+round(size.i/2))
					data.i <- data.i[data.i > 0 & data.i <= L]
					output[i] <- mean(data[data.i])		
				} else {
					data.i <- data.i[data.i > 0 & data.i <= L]
					output[i] <- mean(data[data.i])
				}
		}
	} else {
		for(i in 1:L) {
			data.i <- (i-floor(size/2)):((i-1)+round(size/2))
			data.i <- data.i[data.i > 0 & data.i <= L]
				if(length(data.i) != size) {
					output[i] <- data[i]	
				} else {
					output[i] <- mean(data[data.i])
				}
		}
	}
	return(output)
}

