

source("../volcano_curve_plot.r")
source("../volcano_curve_optim.r")
source("../plot_volcano.r")
source("../plot_volcano_v2.r")
source("../normalise.r")

# really bad fake data
        n = 1000
                data.labels <- sapply(1:n, function(x) {paste(LETTERS[round(runif(5, 1,26))],collapse="")})
                l2fc <- rnorm(n,0,3)
                l2fc <- l2fc[order(abs(l2fc))]
                pval <- sort(-log10(runif(n,0,1)))+rnorm(n,0.5,0.2)
                pval[pval < 0] <- 1e-23
                pval <- 10^-pval
                basemeans <- abs(rnorm(n,500,300))

# data set 
      t1 <- data.frame(data.labels, 
			l2fc,
			pval,
			basemeans
			)

# colours
        colour1 <- runif(3,0,255)
        colour2 <- runif(3,0,255)


# This is a little messy and it isn't one function yet... (It can be compiled as such I guess (see last example))

# First, the volcano curve optimisation (this takes a positive set of genes / features and runs a TPR rate over it.
        vcor <- volcano.curve.optim(t1[,2], -log10(t1[,3]), which(t1[,3] < 0.05), report=T)

# Second, we calculate which genes are above this optimal curve (given by the curve function without plotting)
        vc.signif <- volcano.curve(t1[,2], -log10(t1[,3]), vcor$best.c, colour1, colour2, plot.it=F)

# We then generate the first plot 
        plot.volcano(t1, legend=F, colour1=colour1, colour2=colour2, fade_style=3, p.level=0.00, lines=F)

# We then use the rows from our optimal curve (vc.signif in this case) and plot them as a new layer.
        plot.volcano(t1[vc.signif,], legend=F, colour1=colour1, colour2=colour2, fade_style=3, p.level=1, fold.change=0,lines=F, new.plot=F)

# Finally, I add the curve in using the volcano.curve function and the optimal C
        volcano.curve(t1[,2], -log10(t1[,3]), vcor$best.c, colour1, colour2, plot.it=T, cex=0)




# For example, this combined function could work.

        combined.volcano <- function(dataset, colour1, colour2,...){
        vcor <- volcano.curve.optim(dataset[,2], -log10(dataset[,3]), which(dataset[,3] < 0.05), report=T)
        vc.signif <- volcano.curve(dataset[,2], -log10(dataset[,3]), vcor$best.c, colour1, colour2, plot.it=F)
        plot.volcano(dataset, legend=F, colour1=colour1, colour2=colour2, p.level=0.00, lines=F,...)
        plot.volcano(dataset[vc.signif,], legend=F, colour1=colour1, colour2=colour2, p.level=1, fold.change=0,lines=F, new.plot=F,...)
        volcano.curve(dataset[,2], -log10(dataset[,3]), vcor$best.c, colour1, colour2, cex=0)
        }

combined.volcano(t1, colour1, colour2, fade_style = 3)

# Extra points of interest should be added manually.

# This is not working that well on scales...I have a re-work planned. (see plot volcano v2)
















