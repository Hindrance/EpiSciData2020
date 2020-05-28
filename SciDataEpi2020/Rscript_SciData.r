# The R script for Cytokine Library data analysis 2020.

# Quickly define the platform...
if(length(grep("linux", version)) > 0){
	platform = "linux"
	} else {
	platform = "windows"
}

# Version
analysis.version = "SciData"

# Define directory for output
main.dir = paste("output/v_", analysis.version, sep="")

# meta packages
  library(rmarkdown)
  library(knitr)
  opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

# get new packages?
  GetPackages=F

# source for package loading and personal functions
  source("scripts/packages.r")

######################################################
  version.source("scripts/import_bulk7.r", analysis.version)
######################################################


######################################################
  version.source("scripts/generate_kallisto_estimates.r", analysis.version)
######################################################


######################################################
  version.source("scripts/import_estimates.r", analysis.version)
######################################################


######################################################
  version.source("scripts/technical_validation.r", analysis.version)
######################################################


######################################################
  version.source("scripts/data_processing.r", analysis.version)
######################################################


######################################################
  version.source("scripts/differential_expression_analysis.r", analysis.version)
######################################################


######################################################
  version.source("scripts/interactome_processing.r", analysis.version)
######################################################




# GENERATE FIGURES
######################################################
  version.source("scripts/Figure 2.r", analysis.version)
######################################################


######################################################
  version.source("scripts/Figure 3.r", analysis.version)
######################################################


# Render the report (This report has EVAL=F arguments for a few sections... This means that 
# it compiles much faster as it doesn't need to re-run the analyses code. It uses the objects generated 
# in the scripts above
  rmarkdown::render("report_SciData.Rmd")

