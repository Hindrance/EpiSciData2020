
# Similar function
	similar <- function(data1, data2, round.places=9){
		if(class(data1) == "character"){
			length(which(data1 == data2)) / length(data1) *100

		} else{

		r.v <- round(data1, round.places) == round(data2, round.places)
			similarity.fraction <- length(which(r.v) == T) / length(r.v) * 100
				similarity.fraction
		}
	}
