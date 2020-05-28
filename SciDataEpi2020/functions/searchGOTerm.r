
######################################################
# retrieves information from the EBI resource QuickGO.. 
# uses the API with instructions from:
# "https://www.ebi.ac.uk/QuickGO/api/index.html#/"
# by Vince :D
######################################################

# requires:
  library(httr)
  library(jsonlite)
  library(xml2)

searchGOTerm <- function(keyword, number = 10){
  # requires:
    require(httr)
    require(jsonlite)
    require(xml2)

# Query / search GO for a keyword
    keyword = gsub(" ", "%20", keyword)
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
    for(i in 1:length(downloaded.dataframe)){
      downloaded.dataframe[,i] <- unlist(downloaded.dataframe[[i]])
    }
    colnames(downloaded.dataframe) <- names(fromJSON(json)[[2]])


  return(downloaded.dataframe)
}
  


