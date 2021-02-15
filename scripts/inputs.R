################
#User specified inputs for entire analysis to run
###############

#Implementation settings-------------------------
cores <- 5 #num cores to run parallel processes on (depends on HPC session setup)

#Constants-------------------------------
g <- 9.8 #m/s2

#For BIKER validation---------------------
err <- 1 #set to 1 to corrupt S and dA with Durand et al. layover model
mannings_uncertainity <- 0.25 #see monte_carlo_analysis.R for actual uncertainty
  #For these tests, should assume only error is from manning's equation (0.21), NOT k600 model b/c we are assuming this model is 'true'
  #In actual implementation, we want posterior uncertainity to reflect both (which is ~0.30 per our Monte Carlo Analysis)


#for FCO2 analysis-------------------
molarMass <- 44.01 #g/mol for CO2
pCO2_a <- 390 #uatm
