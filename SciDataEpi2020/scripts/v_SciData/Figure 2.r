
######################################################
# FIGURE 2 - TECHNICAL VALIDATION
######################################################  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# figure A 
 
pdf("Figure 2.pdf", width=10, height=10)
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
 
 
 
