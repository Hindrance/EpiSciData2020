sample.equal <- function(vector, classes, n){
  n2sub <-  min(table(classes))
  if(n < n2sub){a = n} else {a = n2sub}
  
  sub.vector <- unlist(lapply(levels(classes), function(k) {
    sample(which(classes==k), a, replace=F)
    }))
    
  equal.sample.vector <- vector[sub.vector]
  
  return(equal.sample.vector)  
}


