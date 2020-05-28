
# this function is for plotting into single plot figure space. 

subplot.coords <- function(n.y, n.x, bound=c(0,1,0,1)){
    # calculate the square properties of such a plot (rectangular to nearest square OR custom...)
    #    n.x <- floor(sqrt(L))
    #    n.y <- ceiling(sqrt(L))
       
    #sub figure fun2ction (multipanel inside an R plot window...)
    # define coordinates for the subfigure
        x1 = bound[1]
        x2 = bound[2]
        y1 = bound[3]
        y2 = bound[4]
    # Calcula the coordinates of each figure (in reading direction)
        x1.vector <- rep(c((x1+(1:(n.x)-1)*(x2-x1)/n.x)), n.y)
        x2.vector <- rep(c((x1+(1:(n.x)-1)*(x2-x1)/n.x)), n.y)+(x2-x1)/n.x
        y1.vector <- rep(c((y1+((n.y):1-1)*(y2-y1)/n.y)), each=n.x)
        y2.vector <- rep(c((y1+((n.y):1-1)*(y2-y1)/n.y)), each=n.x)+(y2-y1)/n.y
        
        
    # generate a data frame object of our lovely coordinates (we call them by row e.g: coords.df[i,] == c(x1, x2, y1, y2))
        coords.df <- data.frame(x1.vector,x2.vector,y1.vector,y2.vector)
    return(coords.df)
}


