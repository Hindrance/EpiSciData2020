######################################################
######################################################
######################################################




######################################################
# Script details

# NAME
  script.name = "technical_validation"

# DESCRIPTION
  # This is a script with technical validation stuff in. 
  # It has the necessary code and calculations to 
  # generate figure 2 and associated statistics.

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir, "figures")
######################################################



                                                
######################################################
######################################################
######################################################

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# 1  
# Gene / read depth saturation curves / analysis

# Utilise a hypergeometric distribution to subsample our data. 
  library(extraDistr) 
  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# EPI dataset
  # The matrix we want to analyse
  EPI.counts <- EPI.counts[!rownames(EPI.counts) == "null",]
  counts.data.matrix <- as.matrix(EPI.counts)
  n2sample <- round(seq.log.basic(10,max(colSums(counts.data.matrix)), 40))
  set.seed(666)

  res.epi = lapply(1:ncol(counts.data.matrix), function(n){
    z = data.frame(n2sample)
    for(i in 1:length(n2sample)){

    # the number of reads to sample
      k <- n2sample[i]
    
    # the number of samples to create for each sample
      m = 1
    
    # Generate our sub-sampled matrix
      x <- as.numeric(rmvhyper(m,t(as.matrix(counts.data.matrix[,n])),min(k,sum(counts.data.matrix[,n]))))
      
      y  = x > 0
     
      z[i,1:2] = c(min(k,sum(counts.data.matrix[,n])),sum(y))
    }
    z
  })

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# CM dataset  
  # The matrix we want to analyse
  counts.data.matrix <- as.matrix(CM.counts)
  
  n2sample <- round(seq.log.basic(10,max(colSums(counts.data.matrix)), 40))
  set.seed(666)

  res.cm = lapply(1:ncol(counts.data.matrix), function(n){
    z = data.frame(n2sample)
    for(i in 1:length(n2sample)){

    # the number of reads to sample
      k <- n2sample[i]
    
    # the number of samples to create for each sample
      m = 1
    
    # Generate our sub-sampled matrix
      x <- as.numeric(rmvhyper(m,t(as.matrix(counts.data.matrix[,n])),min(k,sum(counts.data.matrix[,n]))))
      
      y  = x > 0
     
      z[i,1:2] = c(min(k,sum(counts.data.matrix[,n])),sum(y))
    }
    z
  })


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# here we make a combined figure for all samples
  gene.y <- sapply(1:length(res.epi), function(i) {res.epi[[i]][,2]})
  depth.x <- sapply(1:length(res.epi), function(i) {res.epi[[i]][,1]})
  x = depth.x
  y = gene.y
  
  plot(as.matrix(x), as.matrix(y), col=rep(c(Discrete[6], Discrete[6], Discrete[6], Discrete[9], Discrete[9], Discrete[9]), each=40),
  xlab="sampled aligned reads (simulated depth)",
#  ylab="Genes counted (with at  5 reads) ", cex=1,
  ylab=expression(paste("Genes counted (with ">=" 1 reads)")), cex=1,
  ylim=c(0,30000), xlim=c(0, 3.2e7), pch=5)
    lines(x[,1], y[,1], col=Discrete[6])
    lines(x[,2], y[,2], col=Discrete[6], lty=2)
    lines(x[,3], y[,3], col=Discrete[6], lty=3)
    
    lines(x[,4], y[,4], col=Discrete[9])
    lines(x[,5], y[,5], col=Discrete[9], lty=2)
    lines(x[,6], y[,6], col=Discrete[9], lty=3)
  text(tail(x[,1],1)*1.1, tail(y[,1],1)+800, c("EPI"), cex=1, col=Discrete[6])
  text(tail(x[,5],1), tail(y[,5],1)-1000, c("NC"), cex=1, col=Discrete[9])
  
  gene.y <- sapply(1:length(res.cm), function(i) {res.cm[[i]][,2]})
  depth.x <- sapply(1:length(res.cm), function(i) {res.cm[[i]][,1]})
  x = depth.x
  y = gene.y
  points(as.matrix(x), as.matrix(y), col=rep(c(Discrete[9], Discrete[9], Discrete[7], Discrete[7]),each=40),
  pch=1, cex=1)
 
  lines(x[,1], y[,1], col=Discrete[9])
  lines(x[,2], y[,2], col=Discrete[9], lty=2)
  lines(x[,3], y[,3], col=Discrete[7])
  lines(x[,4], y[,4], col=Discrete[7], lty=2)
  text(tail(x[,1],1), tail(y[,1],1)-1000, c("D0"), cex=1, col=Discrete[9])
  text(tail(x[,4],1), tail(y[,4],1)-1000, c("D30"), cex=1, col=Discrete[7])

  abline(v = 3e7, lty=2)
  abline(v = 2e7, lty=3)
  
  mtext("A", 3, adj = 0)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# 2  
# LOOKING AT KEY MARKERS OF CELL TYPES 

# Key markers
  CM.markers <- c(
    "MYH6",
    "MYH7",
    "MYL2",
    "MYL7"
  )

  EPI.markers <- c(
    "WT1",
    "TCF21",
    "BNC1",
    "ALDH1A2"
  )

  NC.markers <- c(
    "TFAP2A",
    "SOX9"
   )
 
  H9.markers <- c(
    "POU5F1",
    "SOX2"
  )
  
 
# EPI DATASET   
# Create a very quick size factor calculation...
# Counts per million
  epi.sf <- colSums(EPI.counts+1)/1e6

# and normalise...
  EPI.genes <- t(t(EPI.counts+1)/epi.sf)
  eml <- length(EPI.markers)
  
  EPI.genes = EPI.genes[EPI.markers,]
  
  EPI.genes.df <- data.frame(
      "log.cpm"=log(as.vector(EPI.genes)+1),
      "gene"=rep(rownames(EPI.genes), length(EPI.genes[1,])), 
      "cell.type"=rep(c("EPI", "EPI", "EPI", "NC", "NC", "NC"), each = length(EPI.genes[,1])),
      "group"=paste(rep(rownames(EPI.genes), length(EPI.genes[1,])),rep(c("EPI", "EPI", "EPI", "NC", "NC", "NC"), each = length(EPI.genes[,1])),sep=".")
    )
  EPI.genes.df <- EPI.genes.df[order(EPI.genes.df$group),]

# NC gene expression
  NC.genes <- t(t(EPI.counts+1)/epi.sf)
  nml <- length(NC.markers)
  
  NC.genes = NC.genes[NC.markers,]
  
  NC.genes.df <- data.frame(
      "log.cpm"=log(as.vector(NC.genes)+1),
      "gene"=rep(rownames(NC.genes), length(NC.genes[1,])), 
      "cell.type"=rep(c("EPI", "EPI", "EPI", "NC", "NC", "NC"), each = length(NC.genes[,1])),
      "group"=paste(rep(rownames(NC.genes), length(NC.genes[1,])),rep(c("EPI", "EPI", "EPI", "NC", "NC", "NC"), each = length(NC.genes[,1])),sep=".")
    )
  NC.genes.df <- NC.genes.df[order(NC.genes.df$group),]
  
# combine them:
  EPI.genes.df <- rbind(EPI.genes.df, NC.genes.df)
 
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  
 
# The Kallisto datasets:
  files <- list.files("data/kallisto_out", include.dirs=F, recursive=T, full=T, pattern=".tsv$")[grep("trimmed", list.files("data/kallisto_out", include.dirs=F, recursive=T, full=T, pattern=".h5$"))]
  
# rename the files:
  names(files) <- gsub("(^.*/.*out/|/.*$)", "", files)
# Import the TPM for h9 d0 and d30
  CM.TPM <- sapply(files, function(file) {read.csv(file, stringsAsFactors=F, sep="\t")[,5]})

# adding rownames
  CM.TPM.rownames <- read.csv(files[1], stringsAsFactors=F, sep="\t")[,1]
# cleaning of the names using the previously generated map file tx2gene  
  CM.TPM <- CM.TPM[which(CM.TPM.rownames %in% tx2gene[,1]),]
  CM.TPM.rownames <- CM.TPM.rownames[which(CM.TPM.rownames %in% tx2gene[,1])]
  rownames(CM.TPM) <- tx2gene[match(CM.TPM.rownames, tx2gene[,1]),3]
  
# Cleaning column names
  colnames(CM.TPM) <- as.character(gsub("(^.*/.*out/|/.*$)", "", files))
  
## We will combine the transcript level TPM into gene level TPM by summing the 
## transcripts belonging to each gene...
## make a nice new dataframe for ease
#  CM.counts <- txi.kallisto$counts
#  
#  genelist <-  unique(tx2gene[,3])
#
#  pb <- txtProgressBar(min=0 ,max=length(genelist), style=3)
#  CM.TPM.gene <- t(sapply(1:length(genelist), function(i){
#    setTxtProgressBar(pb, i)
#    return(colSums(CM.TPM[tx2gene[,3] == genelist[i],,drop=F]))
#  }))
#  rownames(CM.TPM.gene) = genelist

  
pb <- txtProgressBar(min=0 ,max=length(CM.markers), style=3)
  CM.TPM.gene <- t(sapply(1:length(c(CM.markers,H9.markers)), function(i){
    setTxtProgressBar(pb, i)
    return(colSums(CM.TPM[rownames(CM.TPM) == c(CM.markers,H9.markers)[i],,drop=F]))
  }))
  rownames(CM.TPM.gene) = c(CM.markers,H9.markers)
  
 
  CM.genes = CM.TPM.gene[CM.markers,]
  cml <- length(CM.markers)
  
  CM.genes.df <- data.frame(
      "tpm"=log(as.vector(CM.genes)+1),
      "gene"=rep(rownames(CM.genes), length(CM.genes[1,])), 
      "cell.type"=rep(c("H9", "H9", "CM", "CM"), each = length(CM.genes[,1])),
      "group"=paste(rep(rownames(CM.genes), length(CM.genes[1,])),rep(c("H9", "H9", "CM", "CM"), each = length(CM.genes[,1])),sep=".")
    )
  
  CM.genes.df <- CM.genes.df[order(CM.genes.df$group),]

### H9 gene expression 
  H9.genes = CM.TPM.gene[H9.markers,]
  hml <- length(H9.markers)
  
  H9.genes.df <- data.frame(
      "tpm"=log(as.vector(H9.genes)+1),
      "gene"=rep(rownames(H9.genes), length(H9.genes[1,])), 
      "cell.type"=rep(c("H9", "H9", "CM", "CM"), each = length(H9.genes[,1])),
      "group"=paste(rep(rownames(H9.genes), length(H9.genes[1,])),rep(c("H9", "H9", "CM", "CM"), each = length(H9.genes[,1])),sep=".")
    )
  
  H9.genes.df <- H9.genes.df[order(H9.genes.df$group),]
  
# Combine them:  
  CM.genes.df <- rbind(CM.genes.df, H9.genes.df)
  
  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##  BAR CHARTS 
# Create our barplot dataframes...
# This is pretty arduous, but it works
# We're just formatting the data tables
  CM.tpm.df <- data.frame(row.names=unique(CM.genes.df[,2]),
      "CM.tpm" = sapply(1:(length(CM.genes.df[,1])/4), function(i){
      a <- seq(1,length(CM.genes.df[,1])-1, by=4)[i]
      b <- seq(2,length(CM.genes.df[,1]), by=4)[i]
      return(mean(CM.genes.df[a:b,1]))
    }),
      "H9.tpm" = sapply(1:(length(CM.genes.df[,1])/4), function(i){
      a <- seq(3,length(CM.genes.df[,1])-1, by=4)[i]
      b <- seq(4,length(CM.genes.df[,1]), by=4)[i]
      return(mean(CM.genes.df[a:b,1]))
    })
  )
  # And one for error bars..?
  CM.error.df <- data.frame(row.names=unique(CM.genes.df[,2]),
      "CM.tpm" = sapply(1:(length(CM.genes.df[,1])/4), function(i){
      a <- seq(1,length(CM.genes.df[,1])-1, by=4)[i]
      b <- seq(2,length(CM.genes.df[,1]), by=4)[i]
      return(sd(CM.genes.df[a:b,1]))
    }),
      "H9.tpm" = sapply(1:(length(CM.genes.df[,1])/4), function(i){
      a <- seq(3,length(CM.genes.df[,1])-1, by=4)[i]
      b <- seq(4,length(CM.genes.df[,1]), by=4)[i]
      return(sd(CM.genes.df[a:b,1]))
    })
  )

# Create our barplot dataframes...
  EPI.tpm.df <- data.frame(row.names=unique(EPI.genes.df[,2]),
      "EPI.tpm" = sapply(1:(length(EPI.genes.df[,1])/6), function(i){
      a <- seq(1,length(EPI.genes.df[,1])-1, by=6)[i]
      b <- seq(3,length(EPI.genes.df[,1]), by=6)[i]
      return(mean(EPI.genes.df[a:b,1]))
    }),
      "NC.tpm" = sapply(1:(length(EPI.genes.df[,1])/6), function(i){
      a <- seq(4,length(EPI.genes.df[,1])-1, by=6)[i]
      b <- seq(6,length(EPI.genes.df[,1]), by=6)[i]
      return(mean(EPI.genes.df[a:b,1]))
    })
  )
  # And one for error bars..?
  EPI.error.df <- data.frame(row.names=unique(EPI.genes.df[,2]),
      "EPI.tpm" = sapply(1:(length(EPI.genes.df[,1])/6), function(i){
      a <- seq(1,length(EPI.genes.df[,1])-1, by=6)[i]
      b <- seq(3,length(EPI.genes.df[,1]), by=6)[i]
      return(sd(EPI.genes.df[a:b,1]))
    }),
      "NC.tpm" = sapply(1:(length(EPI.genes.df[,1])/6), function(i){
      a <- seq(4,length(EPI.genes.df[,1])-1, by=6)[i]
      b <- seq(6,length(EPI.genes.df[,1]), by=6)[i]
      return(sd(EPI.genes.df[a:b,1]))
    })
  )
      
 par(mfrow=c(1,2), cex=0.7)
 # using my custom bar plot function
  bars(EPI.tpm.df, error.data = EPI.error.df, col=c(Discrete[6], Discrete[9]), bar.width=0.31, ylab="log CPM")
    legend("topright", legend=c("EPI", "NC"), col=c(Discrete[6], Discrete[9]), pch=15, bty="n")
    mtext("B", 3, adj=0)
    
  bars(CM.tpm.df, error.data = CM.error.df, col=c(Discrete[7], Discrete[9]), bar.width=0.31, ylab="log TPM")
    legend("topright", legend=c("CM", "H9"), col=c(Discrete[7], Discrete[9]), pch=15, bty="n")  
  
  
# Welches two sample t tests?
  CM.ttest.list <- lapply(unique(CM.genes.df[,2]), function(gene){
    tdf <- CM.genes.df[CM.genes.df[,2] == gene,]
    t.test(tdf[tdf[,3] == "CM",1], tdf[tdf[,3] == "H9",1], paired=F)
  })
  CM.ttests.df = data.frame( unique(CM.genes.df[,2]), Pvalue= sapply(unique(CM.genes.df[,2]), function(gene){
    tdf <- CM.genes.df[CM.genes.df[,2] == gene,]
    unlist(t.test(tdf[tdf[,3] == "CM",1], tdf[tdf[,3] == "H9",1], paired=F)[[3]])
  }))
   
  EPI.ttest.list <- lapply(unique(EPI.genes.df[,2]), function(gene){
    tdf <- EPI.genes.df[EPI.genes.df[,2] == gene,]
    t.test(tdf[tdf[,3] == "EPI",1], tdf[tdf[,3] == "NC",1], paired=F)
  }) 
  EPI.ttests.df = data.frame( unique(EPI.genes.df[,2]), Pvalue= sapply(unique(EPI.genes.df[,2]), function(gene){
    tdf <- EPI.genes.df[EPI.genes.df[,2] == gene,]
    unlist(t.test(tdf[tdf[,3] == "EPI",1], tdf[tdf[,3] == "NC",1], paired=F)[[3]])
  }))
  
  write.table(CM.ttests.df, file=file.path(script.dir, "CM_ttests.csv"), sep=",", quote=F, col.names=T, row.names=F)
  write.table(EPI.ttests.df, file=file.path(script.dir, "EPI_ttests.csv"), sep=",", quote=F, col.names=T, row.names=F)
  
  
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Correlations as per Spearmans: > 0.9
# To conform with encode guidelines for isogenic samples
  
  CM.cors <- cor(CM.counts, method="spearman")
  pdf(file.path(script.dir, "figures", "CM_correlations.pdf"))
 par(mar=c(5,5,5,5), omi=c(1,1,1,1), oma=c(0,0,0,0))
    plot(as.data.frame(CM.counts),main="Spearmans Rho", col="indianred2")
    for(i in 1:16){
    par(fig=as.numeric(subplot.coords(4,4, c(0.1, 0.9, 0.15, 0.95))[i,]))
      mtext(round(t(CM.cors)[i],3),2)
    }
  dev.off()
    
    
  EPI.cors <- cor(EPI.counts, method="spearman")
  pdf(file.path(script.dir, "figures", "EPI_correlations.pdf"))
  par(mar=c(5,5,5,5), omi=c(1,1,1,1), oma=c(0,0,0,0))
    plot(as.data.frame(EPI.counts),main="Spearmans Rho", col="steelblue2")
    for(i in 1:36){
    par(fig=as.numeric(subplot.coords(6,6, c(0.1, 0.9, 0.15, 0.95))[i,]))
      mtext(round(t(EPI.cors)[i],3),2)
    }
  dev.off()
  EPI.cors
  
  
  
  
######################################################

######################################################  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# figure A 
 
pdf(file.path(script.dir, "figures", "Figure 2.pdf"), width=10, height=10)
par(fig=c(0,1,0.3,1),new=F, omi=c(0,0,0,0),oma=c(1,1,1,1), mar=c(5,4,2,1))
  gene.y <- sapply(1:length(res.epi), function(i) {res.epi[[i]][,2]})
  depth.x <- sapply(1:length(res.epi), function(i) {res.epi[[i]][,1]})
  x = depth.x
  y = gene.y
  
  plot(as.matrix(x), as.matrix(y), col=rep(c(Discrete[6], Discrete[6], Discrete[6], Discrete[9], Discrete[9], Discrete[9]), each=40),
  xlab="sampled aligned reads (simulated depth)",
#  ylab="Genes counted (with at  5 reads) ", cex=1,
  ylab=expression(paste("Genes counted (with ">=" 1 reads)")), cex=1,
  ylim=c(0,30000), xlim=c(0, 3.2e7), pch=5)
    lines(x[,1], y[,1], col=Discrete[6])
    lines(x[,2], y[,2], col=Discrete[6], lty=2)
    lines(x[,3], y[,3], col=Discrete[6], lty=3)
    
    lines(x[,4], y[,4], col=Discrete[9])
    lines(x[,5], y[,5], col=Discrete[9], lty=2)
    lines(x[,6], y[,6], col=Discrete[9], lty=3)
  text(tail(x[,1],1)*1.1, tail(y[,1],1)+800, c("EPI"), cex=1, col=Discrete[6])
  text(tail(x[,5],1), tail(y[,5],1)-1000, c("NC"), cex=1, col=Discrete[9])
  
  gene.y <- sapply(1:length(res.cm), function(i) {res.cm[[i]][,2]})
  depth.x <- sapply(1:length(res.cm), function(i) {res.cm[[i]][,1]})
  x = depth.x
  y = gene.y
  points(as.matrix(x), as.matrix(y), col=rep(c(Discrete[9], Discrete[9], Discrete[7], Discrete[7]),each=40),
  pch=1, cex=1)
 
  lines(x[,1], y[,1], col=Discrete[9])
  lines(x[,2], y[,2], col=Discrete[9], lty=2)
  lines(x[,3], y[,3], col=Discrete[7])
  lines(x[,4], y[,4], col=Discrete[7], lty=2)
  text(tail(x[,1],1), tail(y[,1],1)-1000, c("D0"), cex=1, col=Discrete[9])
  text(tail(x[,4],1), tail(y[,4],1)-1000, c("D30"), cex=1, col=Discrete[7])

  abline(v = 3e7, lty=2)
  abline(v = 2e7, lty=3)
  
  mtext("a", 3, adj = -0.05, font=2, padj=-1)
  
 ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
 # figure B 
par(fig=c(0,0.37,0,0.3), new=T, mar=c(4,4,2,1))
 # using my custom bar plot function
  bars(EPI.tpm.df[1:4,], error.data = EPI.error.df[1:4,], col=c(Discrete[6], Discrete[9]), bar.width=0.31, ylab="log CPM", las=3, ylim=c(0,6.5))
    legend("topright", legend=c("EPI", "NC"), col=c(Discrete[6], Discrete[9]), pch=15, bty="n", y.intersp = 0.7)
    mtext("b", 3, adj=-0.05, font=2, padj=-1)
    mtext("Epicardial markers", 3, adj=0.5)
    par(fig=c(0.37,0.5,0,0.3), new=T, mar=c(4,0,2,1))
    bars(EPI.tpm.df[5:6,], error.data = EPI.error.df[5:6,], col=c(Discrete[6], Discrete[9]), bar.width=0.31, ylab="log CPM", las=3, yaxt="n", ylim=c(0,6.5))
    mtext("Neural crest\nmarkers", 3, adj=0.5)
    
  par(fig=c(0.5,0.87,0,0.3), new=T, mar=c(4,4,2,1))
    bars(CM.tpm.df[1:4,], error.data = CM.error.df[1:4,], col=c(Discrete[7], Discrete[9]), bar.width=0.31, ylab="log TPM", las=3, ylim=c(0,10))
      legend("topleft", legend=c("CM", "H9"), col=c(Discrete[7], Discrete[9]), pch=15, bty="n", y.intersp = 0.7)  
    mtext("c", 3, adj=-0.05, font=2, padj=-1)
    mtext("Cardiomyocyte markers", 3, adj=0.5)
    
par(fig=c(0.87,1,0,0.3), new=T, mar=c(4,0,2,1))
  bars(CM.tpm.df[5:6,], error.data = CM.error.df[5:6,], col=c(Discrete[7], Discrete[9]), bar.width=0.31, ylab="log TPM", las=3, yaxt="n", ylim=c(0,10))
  mtext("H9 markers", 3, adj=0.5)
  
  
  
dev.off()
 
# end of script
 
 
 
