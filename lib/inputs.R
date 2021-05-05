################
#User specified inputs for entire analysis to run
###############

###############
#Constants
###############
g <- 9.8 #gravitational acceleration [m/s2]
molarMass <- 44.01 #[g/mol] for CO2
pCO2_a <- 390 #[uatm]
mannings_uncertainity <- 0.25 #see MonteCarlo_analysis.R for actual estimated uncertainty. For BIKER validation, should assume only error is from manning's equation (0.25), NOT k600 model b/c we are assuming this model is 'true'. In actual implementation, we want posterior uncertainty to reflect both
cf <- 0.100 #Henderson found empirically that Cf=0.113, if we use the empirical Strickler relation for n (1923- different rivers), we get Cf=0.09. So, we used Cf=0.100. See the notes!
Sg <- 2.65 #specific gravity

###############
#Implementation settings
###############
cores <- 8 #num cores to run parallel processes on (depends on HPC session setup)

###############
#MCMC analysis setup and munges
###############
M <- 5000 # number of sets of hydraulic 'measurements'
munge <- 0 #only set to 1 to grab M hydraulic measurements. Already done so leave at 0
munge2 <- 1 #Set to 1 to run parallelized MC analyses. These really need to be run on an HPC. This is by far the most compute heavy portion of the project
