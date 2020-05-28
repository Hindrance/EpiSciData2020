
######################################################
# FIGURE 3 - VOLCANO PLOTS
######################################################  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
  interactome.dataframe <- EPI.CM.interactome.cut
# DESEQ2 results visualisation:
# Cardio data
  cmdf <- data.frame(
    "genes" = rownames(CM.res),
    "l2fc" = -CM.res[,2],
    "p.adj" = CM.res[,5],
    "BM" = CM.res[,1]
  )
  cmdf[cmdf[,3]==0,3] <- 1e-310
  CM.receptors <- cmdf[,1] %in% interactome.dataframe[,2]
 
# epi data
  epidf <- data.frame(
    "genes" = rownames(EPI.res),
    "l2fc" = EPI.res[,2],
    "p.adj" = EPI.res[,5],
    "BM" = EPI.res[,1]
  )
  epidf[epidf[,3]==0,3] <- 1e-310
  EPI.secreted <- epidf[,1] %in% interactome.dataframe[,1]
  
# arbitrary gap between volcano dataset points:
  gap=10
  
# We then generate the cm plot, make a report of it

        colour3 <- c("white", Discrete[7])
        colour4 <- c("white", Discrete[9]) 
       t1 <- plot.volcano(cmdf, input.format="custom", legend=F, colourset1=colour3, colourset2=colour4, 
          fade_style=3, p.level=0.00, lines=F, new.plot=F, report=T)
       t1[,"x"] = t1[,"x"]+max(epidf[,2])+gap+abs(min(t1[,1]))
          
          
# We then generate the epi plot, make a report of it
        colour1 <- c("white", Discrete[9])
        colour2 <- c("white", Discrete[6]) 
       t2 <- plot.volcano(epidf, input.format="custom", legend=F, colourset1=colour1, colourset2=colour2, 
          fade_style=3, p.level=0.00, lines=F, new.plot=F, report=T)
       t2[,"x"] = t2[,"x"]
          

pdf("Figure_3_volcano_plot.pdf", width=16, height=9)

  par(cex=2, mar=c(4,4,4,4))
  plot(
    c(min(epidf[,2]), max(epidf[,2])+gap+diff(range(cmdf[,2]))),
    c(0,max(-log10(c(cmdf[,3],epidf[,3])))),
    pch=NA,
    xaxt="n",
    ylab="significance (-log10 P)",
    xlab="log 2 fold change"
  )
  

# Association lines:
  
  for(i in 1:length(interactome.dataframe[,1])){
    x0 = epidf[epidf[,1] == interactome.dataframe[i,1],"l2fc"]
    y0 = -log10(epidf[epidf[,1] == interactome.dataframe[i,1],"p.adj"])
    x1 = t1[rownames(t1) == interactome.dataframe[i,2],"x"]
    y1 = t1[rownames(t1) == interactome.dataframe[i,2],"y"]

    c1 = rgb(colorRamp(colors=c("white", "black"))(  y0 / max(-log10(c(cmdf[,3],epidf[,3]))))/255)
    c2 = rgb(colorRamp(colors=c("white", "black"))(  y1 / max(-log10(c(cmdf[,3],epidf[,3]))))/255)
    
    gradient.line(x0,y0,x1,y1, length=100, col=c(c1,c2), lwd=1, lend=3)
  }
         
         
      
# generate points      
        points(t1[,"x"], -log10(t1[,2]), col=t1$cols, pch=16)
# and overlay the interactome membrane factors 
        points(t1[CM.receptors,"x"], -log10(t1[CM.receptors,2]), bg=t1[CM.receptors,"cols"], pch=21)
  top.genes <- match(interactome.dataframe[,2], rownames(t1[CM.receptors,]))[1:10]
  text.overlay(t1[CM.receptors,"x"][top.genes],-log10(t1[CM.receptors,2])[top.genes], rownames(t1)[CM.receptors][top.genes], cex=0.5, void.x = t1[CM.receptors,"x"], void.y=t1[CM.receptors,"y"])
  mtext("Cardiomyocytes", 3, adj=0.7, col=colour3)
  mtext("H9", 3, adj=0.95, col=colour4)
  
  axis(1,at=seq(-15,15,5)+max(epidf[,2])+gap+abs(min(t1[,1])), labels=seq(-15,15,5))


# We then generate the EPI plot    
        points(t2[,"x"], -log10(t2[,2]), col=t2$cols, pch=16)
# and overlay the interactome membrane factors 
        points(t2[EPI.secreted,"x"], -log10(t2[EPI.secreted,2]), bg=t2[EPI.secreted,"cols"], pch=21)
  top.genes <- match(interactome.dataframe[,1], rownames(t2[EPI.secreted,]))[1:10]
  text.overlay(t2[EPI.secreted,"x"][top.genes],-log10(t2[EPI.secreted,2])[top.genes], rownames(t2)[EPI.secreted][top.genes], cex=0.5, void.x = t2[EPI.secreted,"x"], void.y=t2[EPI.secreted,"y"])
  mtext("Epicardial cells", 3, adj=0.3, col=colour2)
  mtext("Neural Crest", 3, adj=0.05, col=colour1)
  
  axis(1,at=seq(-15,15,5), labels=seq(-15,15,5))
dev.off()


# End of script



