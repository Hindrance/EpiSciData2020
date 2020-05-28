######################################################
######################################################
######################################################




######################################################
# Script details

# NAME
  script.name = "data_processing"

# DESCRIPTION
  # some pre-filtering of data, just to make things a little faster...

# DIRECTORY
  script.dir = file.path(main.dir, script.name)

  mkdir(script.dir)

######################################################



                                                
######################################################
######################################################
######################################################

# We remove all genes 0 reads in all samples...
# creating expressed genes vectors for both
  A.express = names(which(rowSums(EPI.counts) > 0 ))

  B.express = names(which(rowSums(CM.counts)  > 0 ))

# Here's the intersect in case we need it...
  express.intersect = A.express[A.express %in% B.express]

  A = EPI.counts[A.express,]
  B = CM.counts[B.express,]
  
# This was a very simple filtering step. We didn't really need to 
# merge or find the intersect of gene lists because we are interested 
# in fundamentally different genes (cognate receptors in B from factors in A).


