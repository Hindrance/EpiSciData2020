## Linear binning of data in R... Outputs a list of 
# [[1]] A class factor vector
# [[2]] The boundaries used to bin classes (useful for legends etc)

  bin.data <- function(data2fraction, n, bin.fractions = F){
    bin.list <- list()
    bin.classes <- factor(length(data2fraction), levels=1:n)
    if(bin.fractions == F){bin.fractions <- seq(min(data2fraction), max(data2fraction), length=n+1)}
      for(i in n:1){
        #bin.list[[n]] <- which(data2fraction > bin.fractions[n] & data2fraction <= bin.fractions[n+1])
        bin.classes[which(data2fraction <= bin.fractions[i+1])] <- as.character(i)
      }
    bin.classes <- as.factor(bin.classes)
    return(list(bin.classes, bin.fractions))
  }
  
