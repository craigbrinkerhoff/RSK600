################
#User specified inputs for entire analysis to run
###############

###############
#Load/install/check in all packages necessary for analysis
###############
#package check
packages = c("tidyverse", "tidyr", "readr", "parallel", "BIKER", "ncdf4", "lubridate", "segmented", "Metrics", "cowplot", "grid", "RColorBrewer", "colorspace", "ggtext", "ggfortify")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
search()
options(scipen = 999)
theme_set(theme_cowplot())

###############
#Implementation settings
###############
cores <- 5 #num cores to run parallel processes on (depends on HPC session setup)

###############
#Constants
###############
g <- 9.8 #gravitational acceleration [m/s2]
molarMass <- 44.01 #[g/mol] for CO2
pCO2_a <- 390 #[uatm]
mannings_uncertainity <- 0.25 #see MonteCarlo_analysis.R for actual estimated uncertainty
  #For BIKER va.idation, should assume only error is from manning's equation (0.25), NOT k600 model b/c we are assuming this model is 'true'
  #In actual implementation, we want posterior uncertainty to reflect both

###############
#MCMC analysis setup and munges
###############
M <- 5000 # number of sets of hydraulic 'measurements'
munge <- 0 #only set to 1 to grab M hydraulic measurements. Already done so leave at 0
munge2 <- 1 #Set to 1 to run parallelized MC analyses. These really need to be run on an HPC. This is by far the most compute heavy portion of the project
