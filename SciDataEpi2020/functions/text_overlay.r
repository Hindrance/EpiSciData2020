######################################################
#  dependency: wordcloud package 
require("wordcloud")

######################################################

#  x = coordinates of points
#  y = coordinates of points
#  text = labels
#  void.x = coordinates of points (by default)
#  void.y = coordinates of points (by default)
#  freedom = size of point to avoid placing labels over
#  cex = size of text label
#  colour = colour of text label
#  tbs = troubleshooting boolean (tbs = T will show the void points used to avoid point|text overlaps)
#  example = Example == T will plot a small example of the function.

######################################################

  text.overlay <- function(x, y, text, void.x, void.y, freedom=cex, cex=1, colour="black", tbs=F, example=F){
   if(example==T){
      x <- runif(100,0,1)
      y <- runif(100,0,1)
      text <- sapply(1:100, function(i) {paste(sample(letters,7,r=T),collapse="")})
      plot(x,y)
    }
#    x[x == Inf] <- max(x[x != Inf])
#    y[y == Inf] <- max(y[y != Inf])
    
    
    if(missing(void.x)) {void.x = x}
    if(missing(void.y)) {void.y = y}

#    n.border <- 40
#    border.void.x <- c(
#    rep(par("usr")[1],n.border),
#    seq(par("usr")[1], par("usr")[2], length=n.border),
#    rep(par("usr")[2],n.border),
#    seq(par("usr")[1], par("usr")[2], length=n.border)
#    )
#    border.void.y <- c(
#    seq(par("usr")[3], par("usr")[4], length=n.border),
#    rep(par("usr")[3],n.border),
#    seq(par("usr")[3], par("usr")[4], length=n.border),
#    rep(par("usr")[4],n.border)
#    )
#    
    void.x <- c(void.x)#, border.void.x)
    void.y <- c(void.y)#, border.void.y)
    
    mask.colour = rgb(0,0,0,0)
    if(tbs==T){mask.colour = "black"}
    colour.mask = c(rep(mask.colour, length(void.x)), rep(colour, length(x)/length(colour)))
    implottable.space = rep(c(freedom,cex), c(length(void.x), length(x)))#*rep(c(0.5,1), c(length(void.x), length(x)))
    #test <- wordlayout(c(void.x, x), c(void.y,y), c(rep("O",length(void.x)), text), cex=implottable.space)
    textplot(c(void.x, x), c(void.y,y),  c(rep("+",length(void.x)), text), new=F, show.lines=F, col=colour.mask, cex=implottable.space,
    xlim=c(par("usr")[1],par("usr")[2]), ylim=c(par("usr")[3],par("usr")[4]))  
  }
  
  

      
      
      
