#############################################################
# PLOTS: POINTS PLUS DENSITIES - by Juliette (and Vince :D)
#############################################################
#                            __
#                               _.-~  )
#                    _..--~~~~,'   ,-/     _
#                 .-'. . . .'   ,-','    ,' )
#               ,'. . . _   ,--~,-'__..-'  ,'
#             ,'. . .  (@)' ---~~~~      ,'
#            /. . . . '~~             ,-'
#           /. . . . .             ,-'
#          ; . . . .  - .        ,'
#         : . . . .       _     /
#        . . . . .          `-.:
#       . . . ./  - .          )
#      .  . . |  _____..---.._/ _____
#~---~~~~----~~~~             ~~
#
# DENSITY PLOT DOLPHIN ######################################

# datax = x values
# datay = y values
# classes = class factors
# colours = ordered vector of colours (as a HEX string)
# a = plot vertical adjust

plot_plus  <- function(datax, datay, classes, colours, a=10,pch=16,cex=1,ylim, xlim, bw=c(diff(range(datax))/20,diff(range(datay)/20)),...){

  if(class(classes)!="factor"){
	  class.vector <- factor(classes, levels=unique(classes))} else {
	  class.vector <- classes
	}
	class.colours <- colours

# Create a plot range - used for all densities
	
	if(missing(xlim)){x.range <- range(datax)+c((-(range(datax)[2]-range(datax)[1])/10),(range(datax)[2]-range(datax)[1])/10)} else {
	x.range <- range(xlim)+c((-(range(datax)[2]-range(datax)[1])/10),(range(datax)[2]-range(datax)[1])/10)
	}

	if(missing(ylim)){y.range <- range(datay)+c((-(range(datay)[2]-range(datay)[1])/10),(range(datay)[2]-range(datay)[1])/10)} else {
	y.range <- range(ylim)+c((-(range(datay)[2]-range(datay)[1])/10),(range(datay)[2]-range(datay)[1])/10)
	}

# Set up plot window
#	par(mar=c(5,5,5,5),xpd=T,cex=cex)
	par(xpd=T,cex=cex)
	plot(datax,datay,pch=NA,
	     yaxs="i",ylim=y.range,
	     xaxs="i",xlim=x.range,...)

# For each class defined in class list.
	for(i in levels(class.vector)){
		# What data belongs to our class i?
			  i.vector <- which(class.vector == i)
		# What colour belongs to our class i?
			  i.colour <- class.colours[which(levels(class.vector) == i)]
		# Let's plot our class i data
		if(pch < 21){points(datax[i.vector],datay[i.vector],col=i.colour,pch=pch)}
		if(pch >= 21){points(datax[i.vector],datay[i.vector],bg=i.colour,pch=pch, col="black")}

		# Let's find the density of our class y along the x-axis
			  x.density <- density(datax[i.vector],from = min(x.range),to = max(x.range), bw=bw[1])
				x.density$x <- c(min(x.density$x), x.density$x, max(x.density$x))
				x.density$y <- c(0, x.density$y, 0)
                                x.density$y <- x.density$y/length(i.vector)
		# Let's find the density of our class y along the y-axis (here we rename x and y because of reasons)
			  y.density <- density(datay[i.vector],from = min(y.range),to = max(y.range), bw=bw[2])
				y.density$x <- c(min(y.density$x), y.density$x, max(y.density$x))
				y.density$y <- c(0, y.density$y, 0)
                                y.density$y <- y.density$y/length(i.vector)
				  names(y.density)[1:2] <- c("y","x")
		# We need a coefficient to make the density appear on our plot scale - this is scaled to the plot window by "a".
			  x.coeff <- ((y.range[2]-y.range[1])/a)/max(x.density$y)
			  y.coeff <- ((x.range[2]-x.range[1])/a)/max(y.density$x)

		# Let's produce some density plots on the x-axis (axis 3) - adding the y axis coefficient.
			  polygon(x.density$x,x.density$y*x.coeff+max(y.range),col=paste(i.colour,"80",sep=""), border=i.colour)
		# Let's produce some density plots on the y-axis (axis 4) - adding the x axis coefficient.
			  polygon(y.density$x*y.coeff+max(x.range),y.density$y,col=paste(i.colour,"80",sep=""), border=i.colour)
	  	  
	  	  
	  # Plot population distribution values??   - to be continued.. .need to find ways to machine learn / parse population between density distib's
	  #	  x.density$y[which.max(x.density$y)]
	}
par(xpd=F)
box()
}


example.data.function <- function(    class.no = 6, n.factors=2,
					class.length.m = 50,
					class.length.sd = 10,
					norm.m.range = c(-50,50),
					norm.sd.range = c(0,30),
					linear.range = c(-50,50)
  					)     {
# Example data
# class.no = number of classes (clusters)
# class.length.m = number of data in each class
# class.length.sd = s deviation of data

# norm.m.range = lower and upper limits of normal distribution means
# norm.sd.range = lower and upper limits of normal distribution SDs
# linear.range = lower and upper limits of linear function

name.vector <- c(LETTERS,letters,0:9)
multi.factor.set <- t(matrix(rnorm(class.length.m*class.no*n.factors, runif(n.factors, 0, 1), runif(n.factors, 0, 0)), n.factors, class.length.m*class.no))
multi.factor.set <- matrix(rnorm(class.length.m*n.factors, runif(n.factors, 0, 1), runif(n.factors, 0, 0)), n.factors, class.length.m)
class.vector <- character()
class.values.x <- numeric()
class.values.y <- numeric()
class.colours <- character()
	for(i in 1:class.no){
		class.name <- paste(name.vector[round(runif(6,1,length(name.vector)))], collapse="")
		class.length <- round(rnorm(1,class.length.m, class.length.sd))
		class.vector <- c(class.vector, rep(class.name, class.length))
		linearx <- seq(runif(1, linear.range[1], linear.range[2]),runif(1, linear.range[1], linear.range[2]), length=class.length)
		lineary <- seq(runif(1, linear.range[1], linear.range[2]),runif(1, linear.range[1], linear.range[2]), length=class.length)
		class.values.x <- c(class.values.x, (rnorm(class.length,runif(1,norm.m.range[1], norm.m.range[2]),runif(1,norm.sd.range[1], norm.sd.range[2])))+linearx)
		class.values.y <- c(class.values.y, (rnorm(class.length,runif(1,norm.m.range[1], norm.m.range[2]),runif(1,norm.sd.range[1], norm.sd.range[2])))+lineary)
		c.temp <- sample(colours(),1)
		c.temp <- sprintf("#%02X%02X%02X",col2rgb(c.temp)[1], col2rgb(c.temp)[2], col2rgb(c.temp)[3])
		class.colours <- c(class.colours, rep(c.temp,class.length))
	}
	class.vector <- as.factor(class.vector)
	output <- data.frame(class.vector, class.values.x, class.values.y, class.colours, stringsAsFactors=F)
	output
}

if(no.example.plx == F){
#################### EXAMPLE DATA #########################

plot.example.data <- function(){
	temp.data <- example.data.function(5, 50, 10, c(-50, 50), c(0, 30), c(-50, 50))
	plot_plus(temp.data$class.values.x, temp.data$class.values.y, classes = temp.data$class.vector,colours=unique(temp.data$class.colours), pch=21, cex=1.2)
}
plot.example.data()

#################### EXAMPLE DATA #########################
}

