# CMYK colour function... (additive colouring)


cmy = function(c, m, y, alpha, maxColorValue=1){
  if(maxColorValue != 1) { c <- c/maxColorValue; m <- m/maxColorValue; y <- y/maxColorValue }
  c <- 1-c; m <- 1-m; y <- 1-y
  hex <- function(v) substring(rgb(v,0,0),2,3)
  if(!missing(alpha)) alpha <- hex(alpha) else alpha <- ''
  paste0('#',hex(c), hex(m), hex(y), alpha)
}
