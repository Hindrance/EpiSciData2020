######################################################
######################################################
######################################################




######################################################
# Script details

# NAME
  script.name = "import_estimates"

# DESCRIPTION
  # importing the kallisto estimates using the vignette from DESeq2 for estimated count data

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir)

######################################################



                                                
######################################################
######################################################
######################################################

# When we import files from kallisto we are interested in turning our estimated 
# transcript counts into gene-level counts. That means we must collect our estimated 
# transcript IDs into single genes. For context, one gene can have several transcripts 
# as exons or isoforms.
  
# We can perform this using three methods. The best method is perhaps to read the transcript FASTA 
# that we used to make the kallisto index.
# 1. FASTA parsing
# 2. GTF parsing from a known genome GTF file...
# 3. Map transcript IDs using biomaRt - the ensembl R tool for annotations. 

# 1.
# We have indexed our list of transcripts using the compressed fasta. Let's find this
# and parse it.
  cat("\r creating or reading tx2gene.csv ... \n")
  if(length(list.files(script.dir, pattern="tx2gene.csv")) == 0){
# Our FASTA location
  FASTA.loc <- "~/Documents/genomes/human_hg38/Homo_sapiens.GRCh38.cdna.all.fa"
# Using a lazy approach to reading this FASTA and parsing by FASTA entry character: ">"
  FASTA <- read.csv(FASTA.loc, sep=">", stringsAsFactors=F, skip=5, head=F)
  
# calculate fields with the transcript ID
  tx.m <- nchar(FASTA[,2])
# subset by those fields with a transcript ID
  FASTA <- FASTA[which(tx.m > 0),2]
  
# REGEX character substitutions across the field gives us the IDs we want  
  tx.id <- gsub(" cdna.*$", "", FASTA)
  gene.id <- gsub("^.*gene:| gene_bio.*$", "", FASTA)
  gene.symbol <- gsub("^.*gene_symbol:| description.*$", "", FASTA)
  
  tx2gene <- data.frame(tx.id, gene.id, gene.symbol, stringsAsFactors=F)
       # Write file 
      write.table(tx2gene, file=file.path(script.dir, "tx2gene.csv"), row.names=F, col.names=F, sep=",", quote=F)
    } else {
      # Read file 
      tx2gene <- read.csv(file.path(script.dir, "tx2gene.csv"), head=F, stringsAsFactors=F)
  } 
# And there we have it! our tx2gene file used for mapping the transcripts counted in kallisto into genes.
 

# The next step is to read our kallisto estimates in using our map csv. This is 
# part of the recommended workflow from DESeq2 vignette.

# kallisto output directory:
  files <- list.files("data/kallisto_out", include.dirs=F, recursive=T, full=T, pattern=".h5$")[grep("trimmed", list.files("data/kallisto_out", include.dirs=F, recursive=T, full=T, pattern=".h5$"),invert=T)]

# rename the files:
  names(files) <- gsub("(^.*/.*out/|/.*$)", "", files)
# Import using tximport (which applies a counts offset to adjust for transcript estimates or something)
  txi.kallisto <- tximport(files, type = "kallisto",
                            countsFromAbundance="lengthScaledTPM",
                            txOut = F,
                            tx2gene = tx2gene[,c(1,3)]
                            )

# make a nice new dataframe for ease
  CM.counts <- txi.kallisto$counts

# For differential expression analysis, most models assume that we are dealing
# with counts data. Thus we should integerise our dataset.
  CM.counts <- round(CM.counts)

# This is the end of this script... hurrah! We have a nicely mapped









