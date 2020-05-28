volcano.curve <- function(x, y, c, colour1="red", colour2="blue", plot.it=T, pch=21,...){
l2fcs <- x
pvalues <- y
if(class(colour1) == "character"){colour1 <- col2rgb(colour1)}
if(class(colour2) == "character"){colour2 <- col2rgb(colour2)}


	# The sd of x
		x0 = sd(l2fcs)
	# Make x coordinates sequence
		x = seq(x0, max(l2fcs)*10, length = 10000)
	# Now we find y
		y = c / (x - x0)
		
		
	# and we can plot this:
		#curve(c / (x-x0),from = x0, to = 100, n= 100000, add=T, lty=2)
		#curve(c / (-x-x0),from = -100, to = -x0, n= 100000, add=T, lty=2)
	
#		lines(x = -x, y = y, lty=2)
#		lines(x = x, y = y, lty=2)

		#curve(c / (x-x0),from = x0, to = 100, n= 100000, add=T)

# First we take all genes with l2fc > the standard deviation.
	sgv1 <- which(l2fcs < -x0 & pvalues > c / (-l2fcs - x0))
	sgv2 <- which(l2fcs >  x0 & pvalues > c / ( l2fcs - x0))

if(plot.it == T){
	colour1 <- rgb(colour1[1]/255,colour1[2]/255,colour1[3]/255)
	colour2 <- rgb(colour2[1]/255,colour2[2]/255,colour2[3]/255)
	lines(x = -x, y = y, lty=2)
	lines(x = x, y = y, lty=2)
	points(l2fcs[sgv1], pvalues[sgv1], pch=pch, bg=colour1, lwd=2,...)
	points(l2fcs[sgv2], pvalues[sgv2], pch=pch, bg=colour2, lwd=2,...)
} else {
	return(c(sgv1, sgv2)[order(-pvalues[c(sgv1,sgv2)])])}
}


volcano.curve.optim1 <- function(x,y,positive.vector.rows){
l2fcs <- x
pvalues <- y
z <- positive.vector.rows
x0 = sd(l2fcs)
c.init <- min(l2fcs) * (min(abs(na.omit(y))) - x0)
c.max <- max(l2fcs) * (max(abs(na.omit(y))) - x0)
c.seq <- seq(c.init, c.max, length=1000)
c.res <- numeric(1000)
for(i in 1:length(c.seq)){
	c = c.seq[i]
	TP <- length(which(y[z] > c / (x[z] - x0) & y[z] > c / (-x[z] - x0)))
	sgv1 <- which(l2fcs[-z] < -x0 & pvalues[-z] > c / (-l2fcs[-z] - x0))
	sgv2 <- which(l2fcs[-z] >  x0 & pvalues[-z] > c / ( l2fcs[-z] - x0))
	TN <- length(c(sgv1, sgv2))
	#length(which(y[-z] > c / (x[-z] - x0) & y[-z] > c / (-x[-z] - x0)))
	c.res[i] <- TP/(TP+TN)
	
}
	c.seq <- c.seq[c.res != 1]
	c.res <- c.res[c.res != 1]
	cbind("c.value" = c.seq[which.max(c.res)], "max.TPR" = c.res[which.max(c.res)])
}



## Redo for the TPR (TP / (TP + FN)) not (TP / (TP + TN))
## Edited for F1-score

#volcano.curve.optim <- function(x,y,positive.vector.rows, report=F){
#l2fcs <- x
#pvalues <- y
#z.all <- 1:length(x)
#z <- positive.vector.rows
#zN <- (1:length(x))[-z]
#x0 = sd(l2fcs)
## rearranged the equation to derive initial values...
#  c.init <- min(na.omit(pvalues)) * (min(abs(na.omit(l2fcs))) - x0)
#  c.max <- max(na.omit(pvalues)) * (max(abs(na.omit(l2fcs))) - x0)
#c.seq <- seq(c.max,c.init, length=1000)
#c.res <- numeric(1000)
#TP.res <- list()
#TN.res <- list()
#FP.res <- list()
#FN.res <- list()
#conf.matrix <- list()
## pb <- txtProgressBar(0,length(c.seq), style=3)
#for(i in 1:length(c.seq)){
#	c = c.seq[i]

#	#sgv1 <- z[which(l2fcs[z] < -x0 & pvalues[z] > c / (-l2fcs[z] - x0))]
#	#sgv2 <- z[which(l2fcs[z] >  x0 & pvalues[z] > c / ( l2fcs[z] - x0))]

#	#TP <- z[which(y[z] > c / (x[z] - x0) & y[z] > c / (-x[z] - x0))]

#	#sgv1 <- zN[which(l2fcs[zN] < -x0 & pvalues[zN] > c / (-l2fcs[zN] - x0))]
#	#sgv2 <- zN[which(l2fcs[zN] >  x0 & pvalues[zN] > c / ( l2fcs[zN] - x0))]

#	#FP <- c(sgv1, sgv2)

#	#sgv1 <- zN[which(l2fcs[zN] < -x0 & pvalues[zN] < c / (-l2fcs[zN] - x0))]
#	#sgv2 <- zN[which(l2fcs[zN] >  x0 | pvalues[zN] < c / ( l2fcs[zN] - x0))]

#	#FN <- z[!(z %in% TP)]

#	#TN <- z.all[-c(TP,FP,FN)]

#	#length(which(y[-z] > c / (x[-z] - x0) & y[-z] > c / (-x[-z] - x0)))
#	#c.res[i] <- (TP+TN)/(length(z)+TN)




#	sgv1 <- which(x < -x0 & y > c / (-x - x0))
#	sgv2 <- which(x >  x0 & y > c / ( x - x0))
#	t.real <- rep(FALSE, length(z.all)); t.real[z] <- TRUE
#	t.predicted <- rep(FALSE, length(z.all)); t.predicted[c(sgv1, sgv2)] <- TRUE
#	conf.matrix[[i]] <- table(factor(t.real, levels=c("TRUE", "FALSE")), factor(t.predicted, levels=c("TRUE", "FALSE")))
#	TP.res[[i]] <- which(t.real & t.predicted)
#	TN.res[[i]] <- which(!(t.real & t.predicted))
#	FP.res[[i]] <- which(t.predicted & !t.real)
#	FN.res[[i]] <- which(t.real & !t.predicted)
#	# positive recall maximising - c.res[i] <- conf.matrix["TRUE", "TRUE"]/(conf.matrix["TRUE", "TRUE"]+ conf.matrix["TRUE", "FALSE"])
#	# accuracy maximising -	c.res[i] <- (conf.matrix["TRUE", "TRUE"]+conf.matrix["FALSE", "FALSE"])/(length(z.all))
#	# positive predictive value (PPV) * absolute true-positive rate (A.TPR)
##		c.res[i] <- ((conf.matrix[[i]]["TRUE", "TRUE"])/(conf.matrix[[i]]["TRUE", "TRUE"] + conf.matrix[[i]]["FALSE", "TRUE"]))*length(TP.res[[i]])
#	# Secondary deciding factor. Use recall as the factor ONLY if the c.res has more than one maximal point.


#        # F-score
#       
#        c.res[i] <- 2*conf.matrix[[i]]["TRUE", "TRUE"] / 
#        (2*conf.matrix[[i]]["TRUE", "TRUE"] + conf.matrix[[i]]["FALSE", "TRUE"] + conf.matrix[[i]]["TRUE", "FALSE"])
#       # setTxtProgressBar(pb,i)
#   }
#     


#	remove.me <- which(c.res == "NaN")
#	c.res[remove.me] <- 0

##	Primary deciding factor: Highest PPV * A.TPR
##		which(c.res == max(c.res))

## 	Secondary deciding factor. Use recall as the factor ONLY if the c.res has more than one maximal point.
##		which(c.res == max(c.res))[which.max(sapply(which(c.res == max(c.res)), function(i) {length(TP.res[[i]])/length(z)}))]



##	cbind("c.value" = c.seq[which.max(c.res)], "max.TPR" = c.res[which.max(c.res)])

#	#c.seq <- c.seq[c.res != 1]
#	#c.res <- c.res[c.res != 1]
##	cbind("c.value" = c.seq[which.max(c.res)], "max.TPR" = c.res[which.max(c.res)])
##	cbind("c.value" = c.seq[max(which(c.res == max(c.res)))], "max.TPR" = c.res[which.max(c.res)])
##	cbind("c.value" = c.seq[which.min(1-c.res)], "max.TPR" = c.res[which.max(c.res)])

#if(report==T){
#	return(list("TP" = TP.res, "TN" = TN.res, "FP" = FP.res, "FN" = FN.res, "F1" = c.res, 
#			"recalls" = sapply(1:1000, function(i) {length(TP.res[[i]])/length(z)}),
#			"best.res" = which(c.res == max(c.res)),
#			"best.res.recalls" = sapply(which(c.res == max(c.res)), function(i) {length(TP.res[[i]])/length(z)}),
#			"best.c" = c.seq[which.max(c.res)],
#			"best.c.recall" = c.seq[which(c.res == max(c.res))[which.max(sapply(which(c.res == max(c.res)), function(i) {length(TP.res[[i]])/length(z)}))]]
#			)
#	)
#} else {
#	return(cbind("c.value" = c.seq[which(c.res == max(c.res))[which.max(sapply(which(c.res == max(c.res)), function(i) {length(TP.res[[i]])/length(z)}))]], "max.PPV.ATPR" = c.res[which(c.res == max(c.res))[which.max(sapply(which(c.res == max(c.res)), function(i) {length(TP.res[[i]])/length(z)}))]]))
#}
#cat("\n")
#}



