library(httr)
library(jsonlite)
library(xml2)

searchGOTerm <- function(keyword, number = 10){
# Query / search GO for a keyword
    keyword = keyword
  # how many returned ids?
    return.limit = number
  
  # retrieve the URL
    requestURL <- paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=",
                       keyword, 
                       "&limit=", 
                       return.limit, 
                       "&page=1", 
                       sep=""
                       )
                    
    r <- GET(requestURL, accept("application/json"))
  # check for warnings
    stop_for_status(r)
    
  # retrieve JSON formatted file stuff(??)
    json <- toJSON(content(r))
  # convert to actual useful format 
    downloaded.dataframe <- fromJSON(json)[[2]]
  # reformat from lists to character strings:
    downloaded.dataframe <- sapply(1:length(downloaded.dataframe), function(i) {unlist(downloaded.dataframe[[i]])})
    colnames(downloaded.dataframe) <- colnames(fromJSON(json)[[2]])
    
  return(downloaded.dataframe)
}
  








#  library(httr)
#  library(jsonlite)
#  library(xml2)


  GO.terms <- c("GO:0005856")
  
  
  GO.URL <- paste("https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?goId=GO%3A0070125")



requestURL <- "https://www.ebi.ac.uk/QuickGO/search/cytoskeleton"
r <- GET(requestURL, accept("text/tsv"))

#  requestURL <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?goId=GO%3A0070125"
#  r <- GET(requestURL, accept("text/tsv"))

#  stop_for_status(r)

#  results <- content(r,as = "text")


# parses content from the httr package, response(content) afaik...
  
  content.parser <- function(content) {
  # grab our content
    a <- strsplit(results, "\n")
    
#    i = 1
#    b <- strsplit(a[[1]][i],split="\t")
  
  # take the headers, make into valid names
    output.head = make.names(unlist(strsplit(a[[1]][1],split="\t")))
  
  # parse each line individually, make a matrix (sapply) and transpose the matrix.  
    output <- t(sapply(2:length(a[[1]]), function(i) {
       unlist(strsplit(a[[1]][i],split="\t"))
      }  
    ))
  
  # format column names
    colnames(output) <- output.head
    
  # return the output! (wooh!!)  
    return(output)
  }
