# This is an arrow overlay function. It adds an overlay to plots where 
# individual points are linked (I.E time-point data). It is currently limited to 
# two point connections (df[1,n][1:2]) lthough it can be extended by explicitly calling the function 
# for the 2:3, 3:4 and (n):(n+1) cases.

  arrow.overlay <- function(dataset, from=1:(length(dataset[,1])/2), to=(length(dataset[,1])/2+1):(length(dataset[,1])), arrow.factor=0.9, angle = 35, length=0.1, ...){
  # dataset in form of df[x,y]
  arrows(
      dataset[to,1]   - arrow.factor * (dataset[to,1] - dataset[from,1]), dataset[to,2]   - arrow.factor * (dataset[to,2] - dataset[from,2]),
      dataset[from,1] + arrow.factor * (dataset[to,1] - dataset[from,1]), dataset[from,2] + arrow.factor * (dataset[to,2] - dataset[from,2]),
      angle=angle, length=length
    )
  }

