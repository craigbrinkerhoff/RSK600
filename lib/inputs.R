################
#User specified inputs for entire analysis to run
###############

###############
#Constants
###############
g <- 9.8 #gravitational acceleration [m/s2]
molarMass <- 12.01 #[g/mol] for C
pCO2_a <- 400 #[uatm]
uncertainity <- 0.30 #see regression model for actual uncertainty

###############
#Implementation settings
###############
cores <- 31 #num cores to run parallel processes on (depends on HPC session setup)
samplingRate <- 14 #sampling frequency used with the pCO2 data from Beaulieu et al. (2012)
