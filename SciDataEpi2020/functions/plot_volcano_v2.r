# PLOT VOLCANO # by Vincent
require(RColorBrewer)
# This function takes a data frame with three columns and plots a volcano plot!

# Submit data in the form of:

# Column 1: data label text
# Column 2: Log fold changes
# Column 3: p values
# Column 4: base means

# OR submit a DESeq results file and specify the argument: input.format="DESeq"

# Note: colour input has to be in the form of RGB vectors. i.e: c(255, 0, 0) for pure red.

plot.volcano <- function(data, input.format="DESeq", colourset1 = c("white", "blue"), colourset2 = c("white", "red"), 
                          custom.significant, fold.change=2, p.level=0.05, basemean.cut, fade_style, lines, plot.legend, 
                          legendPosition, cex, xlab, ylab, new.plot=T, report = F, do_you_want_text, ...) 
                          {


# Input:
# data = dataframe(dataLabels,logFoldChange,Pvalue,BaseMean) (or see below to use a DESeq results object).

# input.format = can be specified to be a "DESeq" format. For ease of use...

# colourset1 = vector of colours e.g c("white, "blue")

# colourset2 = vector of colours e.g c("white, "red")

# fold.change = FC cut off, (non-log) This is log2d in the function.

# p.level = p value cutoff  (non-log) This is -log10d in the function

# basemean.cut = base mean cut off. This is the 1st quartile by default.

# fade_style = (0:4) 0:1 are fold-change gradient. 2:3 include a cross-fade for p.values. 0 and 2 plot significant values as the max colour in gradient.

# lines = (T / F) includes the cut off lines as dashed lines.

# new.plot = (T / F) to plot a new plot OR just the points..

# report = (T / F) to report the data alongisde calculated values for plotting and the colour scheme for all points

# plot.legend = (T / F) to plot a legend (FALSE by default)

# legendPosition = legend position, "topleft" by default. 

if(missing(fade_style)){fade_style = 0}
if(missing(lines)){do_you_want_lines = F} else {do_you_want_lines = lines}
if(missing(plot.legend)){plot.legend = F}
if(missing(legendPosition)){legendPosition = "topleft"}
if(missing(cex)){cex=2}
if(missing(do_you_want_text)){do_you_want_text = F} else {do_you_want_text = do_you_want_text}
if(missing(xlab)){xlab="log2 fold-change"}
if(missing(ylab)){ylab="-log10 p"}

if(input.format == "DESeq") {data = data.frame("gene" = rownames(data), "L2FC" = data[,2], "p.adj" = data[,6], "basemean" = data[,1])}
if(missing(basemean.cut)){basemean.cut = summary(data[,4])[2]}
# Data manipulation
	# data[is.na(data[,3]),3] <- 1   # REPLACE NA PADJ
	data <- data[!is.na(data[,3]),]  # REMOVE NA PADJ
	
	data1 <- data[data[,2] > 0,] 
	# check for positive fold change genes - create dummy otherwise
	if(length(data1[,1]) != 0) {pos = data1} else {pos <- data.frame(0,0,1,0)}
	pos <- pos[which(as.character(pos[,3]) != "NA"),]

	data2 <- data[data[,2] < 0,]
	# check for negative fold change genes - create dummy otherwise
	if(length(data2[,1]) != 0) {neg = data2} else {neg <- data.frame(0,0,1,0)}
	neg <- neg[which(as.character(neg[,3]) != "NA"),]

	
# remove zeros
	pos[pos[,3] == 0,3] <- min(pos[pos[,3] != 0,3])
	neg[neg[,3] == 0,3] <- min(neg[neg[,3] != 0,3])

	pos[,3] <- -log10(pos[,3])
	neg[,3] <- -log10(neg[,3])

output<- data.frame(
                      data,
                      "x"    = numeric(length(data[,1])),
                      "y"    = numeric(length(data[,1])),
                      "z"    = numeric(length(data[,1])),
                      "x_r"  = numeric(length(data[,1])),
                      "x_r2" = numeric(length(data[,1])),
                      "y_r"  = numeric(length(data[,1])),
                      "z_r"  = numeric(length(data[,1])),
                      "cols" = character(length(data[,1])),
                      stringsAsFactors=F)

     
# Let's go!
###########################################################
  # calculate values and colour gradients
  
if(length(data2[,1]) != 0) {
# negative change
	x <- neg[,2]
	y <- neg[,3]
	z <- neg[,4]
	res.vec1 <- which(y > -log10(p.level) & x < -log2(fold.change) & z > basemean.cut)
  if(!missing(custom.significant)){
    res.vec1 <- match(as.character(custom.significant)[data[custom.significant,2] < 0], rownames(neg))
  }
  
 
# text?
	if(do_you_want_text == T){
	txt.norm.x <- rnorm(length(x), x, (max(x)-min(x))/(cex*50))
	txt.norm.y <- rnorm(length(y), y, (max(y)-min(y))/(cex*50))
	if(length(x[res.vec1]) != 0){text(txt.norm.x[res.vec1], txt.norm.y[res.vec1], neg[res.vec1,1])}}

# colour grading
			x_r2 <- x/(min(x)/2)
			x_r2[x_r2 > 0.6] <- 0.6
			x_r <- x/min(x) 
			y_r <- y/max(y)
			z_r <- x_r * y_r
			z_r <- z_r/max(z_r)
			if(fade_style <= 1){r = x_r} else {r = z_r}
      cols = rgb(colorRamp(colourset1)(normalise(r))/255)
      if(fade_style %in% c(0,2)) {cols[res.vec1] = rgb(colorRamp(colourset1)(1)/255)}
     # if(fade_style <= 3) {cols[res.vec1] = rgb(colorRamp(colourset1)(1)/255)}
      
#output / output separately to remove characters
  output[rownames(neg),c("x", "y", "z", "x_r", "x_r2", "y_r", "z_r")] = cbind(x, y, z, x_r, x_r2, y_r, z_r)
  output[rownames(neg), "cols"] = cols
  
}

####

if(length(data1[,1]) != 0) {
# positive change
	x <- pos[,2]
	y <- pos[,3]
	z <- pos[,4]
	res.vec2 <- which(y > -log10(p.level) & x > log2(fold.change) & z > basemean.cut)
  if(!missing(custom.significant)){
    res.vec2 <- match(as.character(custom.significant)[data[custom.significant,2] > 0], rownames(pos))
  }
  
  
# text?
	if(do_you_want_text == T){
	txt.norm.x <- rnorm(length(x), x, (max(x)-min(x))/(cex*50))
	txt.norm.y <- rnorm(length(y), y, (max(y)-min(y))/(cex*50))
	if(length(x[res.vec2]) != 0){text(txt.norm.x[res.vec2], txt.norm.y[res.vec2], neg[res.vec2,1])}}

# colour grading
			x_r2 <- x/(min(x)/2)
			x_r2[x_r2 > 0.6] <- 0.6
			x_r <- x/min(x) 
			y_r <- y/max(y)
			z_r <- x_r * y_r
			z_r <- z_r/max(z_r)
			if(fade_style <= 1){r = x_r} else {r = z_r}
      cols = rgb(colorRamp(colourset2)(normalise(r))/255)
      if(fade_style %in% c(0,2)) {cols[res.vec2] = rgb(colorRamp(colourset2)(1)/255)}
      
#output / output separately to remove characters
  output[rownames(pos),c("x", "y", "z", "x_r", "x_r2", "y_r", "z_r")] = cbind(x, y, z, x_r, x_r2, y_r, z_r)
  output[rownames(pos), "cols"] = cols
  
}
  
# significant features
# OLD res.vec <- as.numeric(c(rownames(neg)[res.vec1[length(res.vec1):1]], rownames(pos)[res.vec2[length(res.vec2):1]]))
# OLD_edited?  res.vec <- as.numeric(c(res.vec1[length(res.vec1):1], res.vec2[length(res.vec2):1]))
  res.vec <- c(rownames(neg)[res.vec1[1:length(res.vec1)]], rownames(pos)[res.vec2[length(res.vec2):1]])



# Make the plot!
######################################################
if(report==F){
# Initialise plot
	par(cex=cex); if(new.plot==T) {
	  plot(c(min(neg[,2])+min(neg[,2])/15,max(pos[,2])+max(pos[,2])/15),c(0,max(neg[,3],pos[,3])+max(neg[,3],pos[,3])/15), pch=NA, col="white", xlab=xlab, ylab=ylab,
	  ...)
	  }

# plot all
  points(output$x, output$y, pch=16, col=output$cols)

# Plot significant
# OLD  points(output$x[res.vec], output$y[res.vec], pch=21, bg=output$cols[res.vec])
  points(output[res.vec,"x"],output[res.vec,"y"], pch=21, bg=output[res.vec,"cols"])
 
# lines?
	if(do_you_want_lines == T){
	lines(c(-10000,10000),rep(-log10(p.level),2), lty=2, lwd=1)
	lines(c(-log2(fold.change),-log2(fold.change)),c(-10000,10000), lty=2, lwd=1)
	lines(c(log2(fold.change),log2(fold.change)),c(-10000,10000), lty=2, lwd=1)
	}

# Legend?

	if(plot.legend == T){
	par(lend=2); legend(legendPosition, legend=c("Reduced expression", "Increased expression"), pt.bg=c(rgb(colorRamp(colourset1)(1)/255), rgb(colorRamp(colourset2)(1)/255)), lty=0, pch=21, lwd=1, bty="n")
	}
	

  } else {

  rownames(output) <- output[,1]
	  output <- output[,-1]
	  
  return(output)}

##################
# function(data, input.format="DESeq", colourset1, colourset2, custom.significant, fold.change=2, 
#          p.level=0.05, basemean.cut, fade_style, lines, plot.legend, legendPosition, cex, 
 #         xlab, ylab, new.plot=T, report = F, ...)
 
# Input:
# data = dataframe(dataLabels,logFoldChange,Pvalue,BaseMean) (or see below to use a DESeq results object).

# input.format = can be specified to be a "DESeq" format. For ease of use...

# colourset1 = vector of colours e.g c("white, "blue")

# colourset2 = vector of colours e.g c("white, "red")

# fold.change = FC cut off, (non-log) This is log2d in the function.

# p.level = p value cutoff  (non-log) This is -log10d in the function

# basemean.cut = base mean cut off. This is the 1st quartile by default.

# custom.significant = a custom vector of significant rows for data

# fade_style = (0:4) 0:1 are fold-change gradient. 2:3 include a cross-fade for p.values. 0 and 2 plot significant values as the max colour in gradient.

# lines = (T / F) includes the cut off lines as dashed lines.

# new.plot = (T / F) to plot a new plot OR just the points..

# report = (T / F) to report the data alongisde calculated values for plotting and the colour scheme for all points

# plot.legend = (T / F) to plot a legend (FALSE by default)

# legendPosition = legend position, "topleft" by default. 
}

############## EXAMPLE ##############################
# Some dodgy distributions here...
	paste(LETTERS[round(runif(5, 1,26))],collapse="")
		data.labels <- sapply(1:1000, function(x) {paste(LETTERS[round(runif(5, 1,26))],collapse="")})
		logfoldchange <- rnorm(1000,0,2)
		pvalues <- abs(rnorm(1000,0.3,0.5))^2
			pvalues[pvalues > 1] <- 1
		basemeans <- abs(rnorm(1000,500,300))

# Create dataframe with the correct column order...
	t1 <- data.frame(data.labels, 
			logfoldchange,
			pvalues,
			basemeans
			)

# Plot multiple versions...
#	par(mfrow=c(2,2), mar=c(5,5,2,2));

# Here is the function:
#	plot.volcano(t1, cex=1.5)
#	plot.volcano(t1, legend=F, cex=1.5)
#	plot.volcano(t1, legend=F, lines=T, cex=1.5, colour1=c(195,65,173), colour2=c(65,220,205), fade=1)
#	plot.volcano(t1, legend=F, lines=T, text=T, cex=1.5, fade=2); par(mfrow=c(1,1))
######################################################

