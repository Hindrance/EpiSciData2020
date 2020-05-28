######################################################
######################################################
######################################################




######################################################
# Script details

# NAME
  script.name = "differential_expression_analysis"

# DESCRIPTION
  # Here is the script for the differential expression 
  # analysis between both data sets

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir)
  mkdir(script.dir,"figures")
######################################################



                                                
######################################################
######################################################
######################################################

# # # # # # # # # # # # # # # # # # # # # # Questions?

# 1. What genes are upregulated in our EPI dataset vs the NC negative control?
# We first grab the EPI dataset
  sampleTable.A <- data.frame(condition = factor(rep(c("EPI", "NC"), each = 3)))
  sampleTable.A[,1] <- relevel(sampleTable.A[,1], "NC")
  rownames(sampleTable.A) <- colnames(A)

    dds.A <- DESeqDataSetFromMatrix(countData = A, colData = sampleTable.A, design = ~condition)
  
  # Run DESeq2, very simple.
    dds.A <- DESeq(dds.A)
  
  # grab the results of course
    res.A <- results(dds.A, contrast=c("condition",  "EPI", "NC"))
  
  # LFC shrinkage
    resLFC.A <- lfcShrink(dds.A, coef="condition_EPI_vs_NC", type="apeglm")
    resLFC.A <- resLFC.A[order(resLFC.A[,5]),]
    
  # Write output
    write.table(resLFC.A, file=file.path(script.dir, "EPIcounts_resLFCA.csv"), sep=",", col.names=T, row.names=T, quote=F)
    
# 2. What genes are upregulated in the downloaded d30 CM data when compared with their d0 data?
# We grab the CM dataset we prepared
  sampleTable.B <- data.frame(condition = factor(rep(c("D0", "D30"), each = 2)))
    relevel(sampleTable.B[,1], "D0")
  rownames(sampleTable.B) <- colnames(B)

    dds.B <- DESeqDataSetFromMatrix(countData = B, colData = sampleTable.B, design = ~condition)
  
  # Run DESeq2, very simple.
    dds.B <- DESeq(dds.B)
  
  # grab the results of course
    res.B <- results(dds.B, contrast=c("condition",  "D30", "D0"))
  
  # LFC shrinkage
    resLFC.B <- lfcShrink(dds.B, coef="condition_D30_vs_D0", type="apeglm")
    resLFC.B <- resLFC.B[order(resLFC.B[,5]),]
    
  # Write output
    write.table(resLFC.B, file=file.path(script.dir, "EPIcounts_resLFCB.csv"), sep=",", col.names=T, row.names=T, quote=F)
    





  
