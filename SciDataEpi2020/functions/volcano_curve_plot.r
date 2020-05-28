

# Function for plotting the curve based from standard deviation of L2FCs and a positive data set

#       x = lf2cs
#       y = -log10 pvalues

#       colour1 and colour2 (LHS,RHS), (as RGB, same as volcano plot... sorry! 
#       I can add checks and stuff to normalise the colours as in other plot functions but I 
#       just didn't have time yet. (TODO).

#       plot.it = a flag to plot (points + lines) on the active plot window. Set to F for a list 
#       of rows which are positive (above the line) for manual plotting.

#       Tips and tricks: Set cex to 0 if you only want lines (if you were plotting in another way 
#       for example).



volcano.curve <- function(x, y, c, colour1, colour2, plot.it=T, pch=21,...){
# input data
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
