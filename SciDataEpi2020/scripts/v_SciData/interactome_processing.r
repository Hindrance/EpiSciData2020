######################################################
######################################################
######################################################

######################################################
# Script details

# NAME
  script.name = "interactome_processing"

# DESCRIPTION
  # processing workflow for making our secretome interactome

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir)

######################################################



                                                
######################################################
######################################################
######################################################
# Dataset loading
  # Acquire the likely secretome from protein atlas
    secretome = read.csv("https://www.proteinatlas.org/search/protein_class%3APredicted+secreted+proteins?format=tsv", head=T, stringsAsFactors=F, sep="\t")
    secretome.genes <- secretome[,1]
#      write.table(secretome.genes, file="data/secretome_genes.txt", sep="\t", quote=F, col.names=F, row.names=F)
  # Acquire the likely membranome from protein atlas
    membranome = read.csv("https://www.proteinatlas.org/search/protein_class%3APredicted+membrane+proteins?format=tsv", head=T, stringsAsFactors=F, sep="\t")
    membranome.genes <- membranome[,1]
#      write.table(membranome.genes, file="data/membranome_genes.txt", sep="\t", quote=F, col.names=F, row.names=F)
#######################################################

# 1. Find genes that are in stringDB

  # We first load up the stringDB database file using the per-channel scores. This may be interesting as 
  # we can filter the list down by interactions only through experimental evidence.
    stringDB2 <- read.csv(stringsAsFactors=F, head=T, "data/stringDB/9606.protein.links.detailed.v11.0.txt", sep=" ")
    
  # The problem with these stringDB files is that they are written as ENSP (ensembl protein IDs).
  # Luckily, we can access the stringDB info file to convert to gene names and map to the rest of our data.
    stringDB.info <- read.csv(stringsAsFactors=F, head=T, "data/stringDB/9606.protein.info.v11.0_edited.txt", sep="\t")

  # We just need to parse and match this file with our list of StringDB proteins.
    sdb.1 <- stringDB2[,1]
    sdb.2 <- stringDB2[,2]
    sdb.map.1 <- match(sdb.1, stringDB.info[,1])
    sdb.map.2 <- match(sdb.2, stringDB.info[,1])
    
    stringDB2[,1] <- stringDB.info[sdb.map.1,2]
    stringDB2[,2] <- stringDB.info[sdb.map.2,2]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

######################################################
# StringDB filtering - retaining only stringDB interactions that were present in other 
# interactomes or interactions databases
# This is commented out, it was very messy and loads of interactions were removed.
# it probably requires the interactome.psi to be loaded too as per the other 
# commented out script above.
#######################################################
#  version.source("scripts/stringDB2_filtering.r", analysis.version)
#######################################################


######################################################
######################################################
######################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Starting to identify which interactiones we want to keep:

  # We want to grab our differentially upregulated genes from EPI and CM, 
  # filtering them by the secreted and membrane-bound proteins respectively.
  # Let us start with stringDB... It is widely used and the scoring seems 
  # reasonable. 
    score.threshold <- 700
    stringDB.experimental <- stringDB2[stringDB2[,"experimental"] > score.threshold,c(1,2,10)]
    stringDB.full <- stringDB2[,c(1,2,7,10)]
  
   
# Identify those interactions from the epicardial upregulated genes.
  # EPI results dataset
    EPI.res <- resLFC.A
  # remove NAs for our analysis...  These are not good genes.
    EPI.res <- EPI.res[!is.na(EPI.res[,5]),]
  # filter off below p value 1e-2
    EPI.res.cut <- EPI.res[EPI.res[,5] < 1e-2,]
  # Order by LFC as dictated in the nature paper supplement Table 2.  
    EPI.res.cut.lfc <- EPI.res.cut[order(-EPI.res.cut[,2]),]
    
  # Find out which upregulated genes are secreted by EPI...
  # Subset by LFC (positive)
    EPI.secretome <- EPI.res.cut.lfc[EPI.res.cut.lfc[,2] > 0,]
  # filter by secreted genes only
    EPI.secretome <- rownames(EPI.secretome)[rownames(EPI.secretome) %in% secretome.genes]
  # Keep this Dataframe? 
    secreted.genes <- EPI.secretome; length(secreted.genes)# n = 379
  # Match secreted genes with the first column of our interactions database
  # This makes the assumption that as these are secreted proteins.
    EPI.secretome <- EPI.secretome[EPI.secretome %in% stringDB.full[,1]]
  # The fraction of membrane genes that are mapped in stringDB
    length(EPI.secretome)/length(secreted.genes)# 99.5 %
    
  # Find all interactions in stringDB where stringDB is in the EPI.secretome
    EPI.secretome.interactome <- stringDB.full[stringDB.full[,1] %in% EPI.secretome,]
      colnames(EPI.secretome.interactome) <-c("secreted", "interactor", "exp.score", "score")
      
# Identify those interactions from the cardiomyocytes upregulated genes       
  # CM results dataset
    CM.res <- resLFC.B
  # remove NAs for our analysis...  These are not good genes.
    CM.res <- CM.res[!is.na(CM.res[,5]),]
  # filter off below p value 1e-2
    CM.res.cut <- CM.res[CM.res[,5] < 1e-2,]
  # Order by LFC as dictated in the nature paper supplement Table 2.  
    CM.res.cut.lfc <- CM.res.cut[order(-CM.res.cut[,2]),]
   
  # Find out which upregulated genes are membrane-bound by CM...
  # Subset by LFC (positive)
    CM.membranome <- CM.res.cut.lfc[CM.res.cut.lfc[,2] > 0,]
  # filter by membrane genes only  
    CM.membranome <- rownames(CM.membranome)[rownames(CM.membranome) %in% membranome.genes]
  # Keep this Dataframe? 
    membrane.genes <- CM.membranome; length(membrane.genes)# n = 1417
  # match membrane genes with the first column of our interactions database
  # This makes the assumption that as these are membrane proteins.
    CM.membranome <- CM.membranome[CM.membranome %in% stringDB.full[,1]]
  # The fraction of membrane genes that are mapped in stringDB
    length(CM.membranome)/length(membrane.genes)# 95.5 %

  # Find all interactions in stringDB where stringDB is in the CM.membranome
    CM.membranome.interactome <- stringDB.full[stringDB.full[,1] %in% CM.membranome,]
      colnames(CM.membranome.interactome) <- c("membrane", "interactor", "exp.score", "score")

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# We can now find the overlap between the two datasets...
# find our potential crosstalk interactome by:
  # EPI[,1] is EPI secretome
  # EPI[,2] is the interactors with EPI secretome

  # CM[,1] is the CM membranome
  # CM[,2] are the interactors for the CM membranome

  # Therefore.. if we find where EPI[,2] is matched in CM[,1], we will find
  # The CM membrane-bound proteins that interact with EPI-secreted proteins.
  EPI.CM.interactome <- EPI.secretome.interactome[which(EPI.secretome.interactome[,2] %in% CM.membranome.interactome[,1]),]
  EPI.CM.interactome <- EPI.CM.interactome[order(EPI.CM.interactome[,1]),]
  EPI.CM.interactome$BaseMean.EPI <- EPI.res[EPI.CM.interactome[,1],1]
  EPI.CM.interactome$L2FC.EPI <- EPI.res[EPI.CM.interactome[,1],2]
  EPI.CM.interactome$BaseMean.CM <- CM.res[EPI.CM.interactome[,2],1]
  EPI.CM.interactome$L2FC.CM <- CM.res[EPI.CM.interactome[,2],2]
  EPI.CM.interactome$L2FC.average <- (EPI.CM.interactome$L2FC.CM + EPI.CM.interactome$L2FC.EPI) / 2 
  EPI.CM.interactome <- EPI.CM.interactome[order(-EPI.CM.interactome$exp.score),]
    write.table(EPI.CM.interactome, file=file.path(script.dir, "EPI_CM_interactome_full_all_information.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome[,c(1,2,4)], file=file.path(script.dir, "EPI_CM_interactome_full.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome[,1], file=file.path(script.dir, "EPI_interactome_full.txt"), sep="\t", col.names=F, row.names=F, quote=F)
    
    
  EPI.CM.interactome.cut <- EPI.CM.interactome[EPI.CM.interactome[,3] > 700,]
  EPI.CM.interactome.cut <- EPI.CM.interactome.cut[order(-EPI.CM.interactome.cut$L2FC.average),]
    write.table(EPI.CM.interactome.cut, file=file.path(script.dir, "EPI_CM_interactome_cut_all_information.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome.cut[,c(1,2,4)], file=file.path(script.dir, "EPI_CM_interactome_cut.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome.cut[,1], file=file.path(script.dir, "EPI_interactome_cut.txt"), sep="\t", col.names=F, row.names=F, quote=F)


  EPI.CM.interactome.ext <- EPI.CM.interactome[which(
    EPI.CM.interactome[,1] %in% as.vector(unlist(EPI.CM.interactome.cut[,1:2])) & 
    EPI.CM.interactome[,2] %in% as.vector(unlist(EPI.CM.interactome.cut[,1:2]))
  ),]
  
  EPI.CM.interactome.ext <- EPI.CM.interactome.ext[EPI.CM.interactome.ext[,3] > 400,]  
    write.table(EPI.CM.interactome.ext, file=file.path(script.dir, "EPI_CM_interactome_ext_all_information.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome.ext[,c(1,2,4)], file=file.path(script.dir, "EPI_CM_interactome_ext.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome.ext[,1], file=file.path(script.dir, "EPI_interactome_ext.txt"), sep="\t", col.names=F, row.names=F, quote=F)

# We are interested to refine this list, perhaps retaining only those genes
# which are secreted by the EPI, and not those which are also secreted by the 
# CMs. To do this we attempt to compare between experimental groups, 
# looking at counts per million read as a normalised between group approach. 
  sec.genes <- unique((EPI.CM.interactome.cut[,1]))
    
  # First we create the CPM matrices
  # CM DATASET   
  # Create a very quick size factor calculation...
  # Counts per million
    cm.sf <- colSums(CM.counts+1)/1e6

  # and normalise...
    CM.cpm <- t(t(CM.counts+1)/cm.sf)
    
    CM.cpm.sec = CM.cpm[sec.genes,3:4]
    

  # EPI DATASET   
  # Create a very quick size factor calculation...
  # Counts per million
    epi.sf <- colSums(EPI.counts+1)/1e6

  # and normalise...
    EPI.cpm <- t(t(EPI.counts+1)/epi.sf)
    
    EPI.cpm.sec = EPI.cpm[sec.genes,1:3]


# Take a log of the CPM (logCPM)
    EPI.logcpm.sec <- log(EPI.cpm.sec + 1)
    CM.logcpm.sec  <- log(CM.cpm.sec  + 1)

# CPM means dataframe 
# Dataframe for results
  cpm.df <- data.frame(row.names=sec.genes, rowMeans(EPI.logcpm.sec), rowMeans(CM.logcpm.sec))
    cpm.order <- order(-(cpm.df[,1]  - cpm.df[,2]))
  cpm.df <- cpm.df[cpm.order,]  

# sd CPM error bars dataframe    
  error.df <- data.frame(row.names=sec.genes, rowSds(EPI.logcpm.sec), rowSds(CM.logcpm.sec))
    error.df <- error.df[cpm.order,]  

# t tests for our datasets
# We can't really be sure of the distributions though... 
# Certainly because the n = 3 and n = 2
  ttest.p <- data.frame(row.names=sec.genes[cpm.order], 
              -(cpm.df[,1]  - cpm.df[,2]),
              round(sapply(1:length(sec.genes), function(i) {
              t.test(EPI.logcpm.sec[i,], CM.logcpm.sec[i,], paired=F)$p.value
              }),3)[cpm.order]
             )
# Threshold    for t test cutoff          
  ttest.cut <- ttest.p[,2] < 0.05

# A list of all genes which are significantly more expressed in EPI
  epi.sec <- rownames(ttest.p)[ttest.p[,2] < 0.05 & ttest.p[,1] < 0]

# subset the interactome results table by the significantly "more EPI" genes
  EPI.CM.interactome_EPI.secreted <- EPI.CM.interactome.cut[EPI.CM.interactome.cut[,1] %in% epi.sec,]
    write.table(EPI.CM.interactome_EPI.secreted, file=file.path(script.dir, "EPI_CM_interactome_EPI_secreted_all_information.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome_EPI.secreted[,c(1,2,4)], file=file.path(script.dir, "EPI_CM_interactome_EPI_secreted.csv"), sep=",", col.names=T, row.names=F, quote=F)
    write.table(EPI.CM.interactome_EPI.secreted[,1], file=file.path(script.dir, "EPI_interactome_EPI_secreted.txt"), sep="\t", col.names=F, row.names=F, quote=F)
  
 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## BAR CHARTS
  
# dataframe of x and y coordinates (results from bars function)     #  
  bar.coords <- bars(cpm.df, error.data = error.df,
    col=c(Discrete[6], Discrete[7]), bar.width=0.31, ylab="log CPM", xaxt="n")
# legend
  legend("topright", legend=c("EPI", "CM"),
    col=c(Discrete[6], Discrete[7]), pch=15, bty="n")

# Axis labels
  axis(1, at=c(1:length(cpm.df[,1])), labels=rownames(cpm.df), las=2)

# Error points
  points(bar.coords[rep(ttest.cut, each=2),1], bar.coords[rep(ttest.cut, each=2),2]+0.6, pch="*")

# line to separate more or less expressedin EPI 
  abline(v=which(ttest.p[,1] > 0)[1]-0.5, lty=2)

#  order:
#  RSPO1     LGR5
#  RSPO1     LGR4
#  RSPO1     RNF43
#  
#  NTF4      NTRK2
#  NTF4      BDNF
#  
#  LY96      TLR4
#  
#  TGFB1     ITGB6
#  TGFB1     ITGAV
#  
#  BMP4      BMPR1B
#  BMP4      ACVR2A
#    
#  IGF2      CLN5
#  
#  INS-IGF2  SPN

#  IL6ST     IL6R
#  
#  APOC1     MPC2
#  
#  PVR       TIGIT 
#  
#  PCDHA10   PCDHA9
#  PCDHA10   PCDHA6 
#  
#  TIMP1     MMP14
#  
#  AGRN      LRP4
#  
#  
#  
#  
#  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
