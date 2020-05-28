seq.log.basic <- function(lower, upper, length) {
		10^seq(log10(lower), log10(upper),length=length)
		#as.numeric(paste(seq(1-(1-by), 1, by=by), "e", rep(seq(lower,upper), each=1/by), sep=""))
		#as.numeric(paste(seq(1, 10-by, by=by), "e", rep(seq(lower,upper), each=(10/by)-(1/by)), sep=""))
	}

