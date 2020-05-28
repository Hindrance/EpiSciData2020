######################################################
######################################################
# I made this as I was fed up of colouring individually by gene each time I was search for expression.

  gene.colour <- function(X, genes, colours=alt.cols, method="mean"){
  # This function is to take an input in the form of gene vector, a gene matrix X[genes, cells], and a colour gradient
  # The output is a list of colours which belong to each cell / sample
  
  # first, we create the gene values (in the case of a 1 gene vector, the sum is the same as the single value...).
  
#    x <- X[genes,,drop=F]
#    
#    for(i in 1:length(genes)){
#    x[i,] <- normalise(x[i,])
#    }
#    
#    for(i in 1:length(genes)){
#    x[i,] <- x[i,] - mean(x[i,])
#    }
#    
    gene.values <- colSums(X[genes, , drop=F])
    if(method=="mean") {gene.values <- colMeans(X[genes, , drop=F])}
  # Then we simply cut it linearly to create a gradient.
    gene.col <- colorRampPalette(color=colours)(100)[as.numeric(cut(gene.values,100))]
    return(gene.col)
  }

