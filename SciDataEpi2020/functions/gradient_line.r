# This function generates a gradient between two points (linear gradient, linear line)  
  
  
  gradient.line <- function(x0, y0, x1, y1, length=10, cols=c("white", "black"), ...) {
    l = length
    x.vector <- sapply(0:l, function(i) {x0+diff(c(x0,x1)) * ((0+i) / l)})
    y.vector <- sapply(0:l, function(i) {y0+diff(c(y0,y1)) * ((0+i) / l)})

    x0 = x.vector[(1:length(x.vector))-1]
    x1 = x.vector[2:length(x.vector)]
    y0 = y.vector[(1:length(y.vector))-1]
    y1 = y.vector[2:length(y.vector)]

    segments(x0,y0,x1,y1, col=colorRampPalette(colors=cols)(l), ...)
  }

  
