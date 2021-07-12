################
#User specified inputs for entire analysis to run
###############

###############
#Constants
###############
g <- 9.8 #gravitational acceleration [m/s2]
molarMass <- 12.01 #[g/mol] for C
pCO2_a <- 390 #[uatm]
uncertainity <- 0.20 #see regression model for actual uncertainity

###############
#Implementation settings
###############
cores <- 8 #num cores to run parallel processes on (depends on HPC session setup)
