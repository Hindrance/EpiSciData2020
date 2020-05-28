
# This is the first function I ever wrote in R. Probably my most used in modelling when taking data from literature or 
# comparing concentrations 

convert <- function(x, x_molar_mass_Da, mass_unit, volume_unit, molarity){
	mass_vector<-c("g", "mg", "ug", "ng", "pg")
	mass_vector_values<-c(1, 1000, 1000000, 1000000000, 1000000000000)
	volume_vector<-c("l", "dl", "ml", "ul")
	volume_vector_values<-c(1, 10, 1000, 1000000)
	molarity_vector<-c("M", "mM", "uM", "nM", "pM")
	molarity_vector_values<-c(1, 1000, 1000000, 1000000000, 1000000000000)
	
	for(i in 1:5){
		if(mass_unit == mass_vector[i]) to_g <- mass_vector_values[i]
	}
	for(i in 1:4){
		if(volume_unit == volume_vector[i]) to_l <- volume_vector_values[i]
	}
	for(i in 1:5){
		if(molarity == molarity_vector[i]) to_desired_conc <- molarity_vector_values[i]
	}
	x <- x/to_g
	x <- x*to_l
	x <- x/x_molar_mass_Da
	x <- x*to_desired_conc
	x
}
