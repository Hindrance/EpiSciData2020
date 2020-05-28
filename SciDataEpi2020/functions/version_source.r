


version.source <- function(default.path, scripts.version, version.control=T,report=T){
		    
	scripts.path = paste("scripts/v_",scripts.version,sep="")
	# errors! Make sure you make that scripts directory :D
	if(!dir.exists(file.path(scripts.path))){
	#cat(paste("\n\rSCRIPTS VERSION ", scripts.version, " NOT READY: NO SCRIPTS DIRECTORY FOUND.\n 
		#Please add a scripts directory labelled as scripts/v_",scripts.version," and fill it with scripts\n", sep=""))
		stop(paste("\n\rSCRIPTS VERSION ", scripts.version, " NOT READY: NO SCRIPTS DIRECTORY FOUND.\n 
		Please add a scripts directory labelled as scripts/v_",scripts.version," and fill it with the scripts\n", sep=""),call.=F)
	} else {
	if(report==T){cat(paste("\rSCRIPTS VERSION ", scripts.version, " READY: SCRIPTS DIRECTORY FOUND AT scripts/v_",scripts.version,"\n\n", sep=""))}
	}
		    
	if(version.control==T){
		new.path = gsub("scripts/", paste(scripts.path,"/",sep=""), default.path)
		source(paste(new.path, sep=""))
	}
	if(version.control==F){
		source(paste(default.path, sep=""))
	}
}

