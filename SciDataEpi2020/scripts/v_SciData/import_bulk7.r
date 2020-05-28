######################################################
######################################################
######################################################




######################################################
# Script details

# NAME
  script.name = "import_bulk7"

# DESCRIPTION
  # Small script to import the previous bulk7 data

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir)

######################################################

cat("\n \r... reading bulk7 data \n")                                              
######################################################
######################################################
######################################################

# Begin generally with importing the data... easy enough
  bulk7 <- read.csv("data/bulk7/bulk7_counts.csv", head=T, row.names=1, stringsAsFactors=F)
  
# grab the epi and NC only columns...
  EPI.counts <- bulk7[,grep("EPI.$|NC", colnames(bulk7))]

# We have gene symbols as the rownames so we will try to map them to each other after...
# Oh well, that was a simple script!


# END
