# returns a GO ID dataframe..

  returnGO <- function(GOIDs) {
    require(httr)
    require(jsonlite)
    require(xml2)

    quickGOquery <- paste(GOIDs, collapse="%2C")
    quickGOquery <- gsub(":", "%3A", quickGOquery)
    quickGOquery <- paste("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", quickGOquery, sep="")
    r <- GET(quickGOquery, accept("application/json"))

    stop_for_status(r)

    json <- toJSON(content(r))
    dl.output <- fromJSON(json)
    GO.list <- data.frame(
      unlist(dl.output$results$id),
      unlist(dl.output$results$name),
      unlist(dl.output$results$definition[,1]),
      stringsAsFactors=F
      )
    return(GO.list)
  }
  
  
  
  
