
## suppose I have a list of counts with rownames that I want to change   
#   bulk.counts <- read.csv("../misc/Bulk_seq_seven/data/bulk-counts.csv", head=T, sep=",", stringsAsFactors=F)
#    row.names(bulk.counts) <- bulk.counts[,1]  
#   
#  names2change <- rownames(bulk.counts)
#  
## We can update all names by applying the following script:

## The URL of our table? Columns for this one:
#map.table <- "https://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&col=gd_pub_refseq_ids&col=md_refseq_id&col=md_ensembl_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"
#    # 1. "Approved.symbol"
#    # 2. "Approved.name"
#    # 3. "Previous.symbols"
#    # 4. "Synonyms"  
#    # 5. "Ensembl.gene.ID"
#    # 6. "RefSeq.IDs"
#    # 7. "RefSeq.supplied.by.NCBI." 
#    # 8. "Ensembl.ID.supplied.by.Ensembl."


## Our mapping table (load that URL and parse)
#      HGNC.map <- read.csv(map.table, sep="\t", head=T, stringsAsFactors=F) 

#      columns2parse <- which(colnames(HGNC.map) %in% c("Approved.name", "Previous.symbols", "Synonyms"))

## Length of list (the number of approved HGNC symbols)
#    L.m <- length(HGNC.map[,1])

## A map matrix - a rapid way of performing matches.
#    HGNC.matrix <- matrix(character(), L.m, 50)
#            
## Set up the mapping matrix
#    pb <- txtProgressBar(min=0 ,max=L.m, style=3)
#    
## For each gene...    
#  for(j in 1:L.m){
#    # generate a vector of unparsed symbols and parsed character strings (for the name, previous names and synonyms)...
#      v.j <- as.character(c(unlist(HGNC.map[j,]),unlist(strsplit(unlist(HGNC.map[j,columns2parse]),split=", | "))))
#    # Update the matrix row with all of these values
#      HGNC.matrix[j,1:length(v.j)] = v.j

#      setTxtProgressBar(pb, j)
#  }

## This gives us a neat representation of all search fields. We can then perform matches to them (where the matrix can 
## be represented by a vector). This is fast because we perform a simple match function across this vector and calculate 
## the matrix row position by calculating the remainded after dividing by the length of rows.
## Pros:
##    - fast
##    - matches FIRST hit (Critical as it will match the approved symbol first).

## Match locations (where are our names first found in this vector?)
#  mid <- match(names2change, HGNC.matrix)

## We divide these positions by the length of matrix + 1 and round DOWN. This leaves a remainder. The reason that we add one to the length is to 
## fix cases where the match position falls exactly on the end of a column. Leaving a remainder of 0.
## (where we know that it belongs to end of the previous column, not the the next column)

## We multiply the rounded factor by the length and subtract this from the match position to achieve the final position of the row. Simple right?
#  mid2 <- mid -  L.m * floor(mid / (L.m+1))
#     
## Now we identify which rows are which, we can map these names back by simply calling the row position and choosing ANY column we want from our map 
## matrix
#  newnames <- HGNC.matrix[mid2,1]

## compare with our old names??? 
#  data.frame(newnames, names2change)






############################################################################################################
############################################################################################################
############################################################################################################


# All in one function?

  update_labels <- function(
                          # The list of names we want to change
                            names2change, 
                          # The features we want to search for our current name against
                            Features = c("Approved symbol", 
                                         "Approved name", 
                                         "Previous symbols", 
                                         "Synonyms",
                                         "Ensembl gene ID",
                                         "RefSeq IDs",
                                         "RefSeq",
                                         "Ensembl ID"
                                         ),
                          # The new feature label that we want to assign
                            newlabels = "Approved symbol",
                          
                          # Progress bar?
                            bar=F
                            )
  {
  # Default:
    # 1. "Approved.symbol"
    # 2. "Approved.name"
    # 3. "Previous.symbols"
    # 4. "Synonyms"  
    # 5. "Ensembl.gene.ID"
    # 6. "RefSeq.IDs"
    # 7. "RefSeq.supplied.by.NCBI." 
    # 8. "Ensembl.ID.supplied.by.Ensembl."

# To make this into a function we need to consider usability and custom requests...
# We can access HGNCs database using custom URLs.

  HGNC.Labels <- data.frame(
    "DataRequest" = c("HGNC ID", "Approved symbol", "Approved name", "Status", "Locus type", "Locus group",
                      "Previous symbols", "Previous name", "Synonyms", "Name synonyms", "Chromosome",
                      "Date approved", "Date modified", "Date symbol changed", "Date name changed",
                      "Accession numbers", "Enzyme IDs", "NCBI Gene ID", "Ensembl gene ID",
                      "Mouse genome database ID", "Specialist database links", "Specialist database IDs",
                     "Pubmed IDs", "RefSeq IDs", "Gene group ID", "Gene group name", "CCDS IDs", "Vega IDs",
                      "Locus specific databases", "NCBI Gene ID", "OMIM ID", "RefSeq", "UniProt ID", "Ensembl ID",
                      "Vega ID", "UCSC ID", "Mouse genome database ID", "Rat genome database ID",
      stringsAsFactors=F
    ),
    "HGNC_label" = c("gd_hgnc_id", "gd_app_sym", "gd_app_name", "gd_status", "gd_locus_type", "gd_locus_group",
                     "gd_prev_sym", "gd_prev_name", "gd_aliases", "gd_name_aliases", "gd_pub_chrom_map",
                     "gd_date2app_or_res", "gd_date_mod", "gd_date_sym_change", "gd_date_name_change",
                     "gd_pub_acc_ids", "gd_enz_ids", "gd_pub_eg_id", "gd_pub_ensembl_id", "gd_mgd_id", "gd_other_ids",
                     "gd_other_ids_list", "gd_pubmed_ids", "gd_pub_refseq_ids", "family.id", "family.name", "gd_ccds_ids",
                     "gd_vega_ids", "gd_lsdb_links", "md_eg_id", "md_mim_id", "md_refseq_id", "md_prot_id", "md_ensembl_id",
                     "md_vega_id", "md_ucsc_id", "md_mgd_id", "md_rgd_id",
      stringsAsFactors=F
    ) 
  )

  request.list <- as.character(HGNC.Labels[match(Features, HGNC.Labels[,1]),2])

cat("\r Requesting dataframe from HGNC server... \n")
# compile the request automatically...
  map.table <- paste("https://www.genenames.org/cgi-bin/download?", 
    paste("col=",request.list, "&", sep="", collapse=""),  
    "status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit", 
    sep=""
  )
  
  
# Our mapping table (load that URL and parse)
      HGNC.map <- read.csv(map.table, sep="\t", head=T, stringsAsFactors=F)
      
      columns2parse <- which(colnames(HGNC.map) %in% c("Approved.name", "Previous.symbols", "Synonyms"))
      newlabel.col <- match(newlabels, Features)
      
# Length of list (the number of approved HGNC symbols)
    L.m <- length(HGNC.map[,1])

# A map matrix - a rapid way of performing matches.
    HGNC.matrix <<- matrix(character(), L.m, 50)

cat("\r Generating mapping matrix.... \n")            
# Set up the mapping matrix
    if(bar==T){pb <- txtProgressBar(min=0 ,max=L.m, style=3)}
    
# For each gene...    
  for(j in 1:L.m){
    # generate a vector of unparsed symbols and parsed character strings (for the name, previous names and synonyms)...
      v.j <- as.character(c(unlist(HGNC.map[j,]),unlist(strsplit(unlist(HGNC.map[j,columns2parse]),split=", | "))))
    # Update the matrix row with all of these values
      HGNC.matrix[j,1:length(v.j)] <- v.j

      if(bar==T){setTxtProgressBar(pb, j)}
  }

# This gives us a neat representation of all search fields. We can then perform matches to them (where the matrix can 
# be represented by a vector). This is fast because we perform a simple match function across this vector and calculate 
# the matrix row position by calculating the remainded after dividing by the length of rows.
# Pros:
#    - fast
#    - matches FIRST hit (Critical as it will match the approved symbol first).

  newnames <- character(length(names2change))
# Match locations (where are our names first found in this vector?)
  mid <- match(names2change, HGNC.matrix)

# We divide these positions by the length of matrix + 1 and round DOWN. This leaves a remainder. The reason that we add one to the length is to 
# fix cases where the match position falls exactly on the end of a column. Leaving a remainder of 0.
# (where we know that it belongs to end of the previous column, not the the next column)

# We multiply the rounded factor by the length and subtract this from the match position to achieve the final position of the row. Simple right?
  mid2 <- mid -  L.m * floor((mid-1) / (L.m))
     
# Now we identify which rows are which, we can map these names back by simply calling the row position and choosing ANY column we want from our map 
# matrix
  newnames[is.na(mid2)] <- names2change[is.na(mid2)]
  newnames[!is.na(mid2)] <- HGNC.matrix[as.numeric(na.omit(mid2)),newlabel.col]
  newnames[newnames==""] <- names2change[newnames==""]
  
# compare with our old names??? 
  return(newnames)
  }
  
  
  
  
