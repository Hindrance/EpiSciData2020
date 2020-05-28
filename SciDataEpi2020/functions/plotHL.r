# This is a plot half-life function. It takes a line and plots the values of the half-life across it... Simple right?
# A very early (poorly written) function of mine.

#########################
#########################
######## plot HL ########
#########################
#########################
plotHL <- function(data, Time_col_num, Conc_col_num, count, scale){
#Calls plot HL with arguments: data (data set object name), Time_col_num (Column number of object that represents time / x), Conc_col_num (column number of object that represents concentration / y), count (number of half-life values to plot - starts from highest conc value).
colnames(data)[Time_col_num] <- "times"
colnames(data)[Conc_col_num] <- "concentrations"
#data[1,Conc_col_num] <- max(data[,Conc_col_num])
if(scale == "log") plot(data[,Time_col_num],data[,Conc_col_num], xlim=c(min(data[,Time_col_num]), max(data[,Time_col_num])), type="l", ylim=c(min(data[,Conc_col_num])+(1e-12), max(data[,Conc_col_num])), main="concentration over time", ylab="Concentration", xlab="Time", log="y")
if(scale == "linear") plot(data[,Time_col_num],data[,Conc_col_num], xlim=c(min(data[,Time_col_num]), max(data[,Time_col_num])), type="l", ylim=c(min(data[,Conc_col_num]), max(data[,Conc_col_num])), main="concentration over time", ylab="Concentration", xlab="Time")
legend("topright", legend=c("Species", "half-life"), lty=c(1, 0), pch=c(NA, "x"), col=c("black", "blue"))
times <<- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
values <- c(max(data[,Conc_col_num]),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
half_lives <- data.frame("Time Interval" = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), "Concentration" = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), "Half Life" = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
#count = 25
for(i in 2:count){
values[i] <- values[i-1]/2
}
for(i in 1:count){
y = as.numeric(mean(subset(data, select=c(concentrations), concentrations > values[i]-(values[i]/10) & concentrations < values[i]+(values[i]/10))[,1]))
x = as.numeric(mean(subset(data, select=c(times), concentrations > values[i]-(values[i]/10) & concentrations < values[i]+(values[i]/10))[,1]))
#lines(c(0.000001, x), c(y,y),col="blue")
#lines(c(x, x), c(y,0.000001),col="blue")
times[i] = x
#text((((times[i+1]+times[i])/2)), 0.001, paste(round(times[i+1]-times[i],3)), col="red", cex = 0.5)
}
for(i in 1:count){
if(times[i+1] != 0) (text(((times[i+1]+times[i])/2), (((values[i]+values[i+1])/2))+(((values[i]+values[i+1])/2))/1, paste(round(times[i+1]-times[i],3)), col="blue", cex = 1))
if(times[i+1] != 0) half_lives[i,1] = times[i]-times[1]
if(times[i+1] != 0) half_lives[i,2] = values[i]
if(times[i+1] != 0) half_lives[i,3] = paste(times[i+1]-times[i])
half_lives <<- half_lives
#half_lives_ <<- assign(paste0("half_lives_", deparse(substitute(data))), half_lives)
}
}


#########################
#########################
#########################
#########################
#########################

