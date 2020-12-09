#Runs the three analysis scripts for this project.
#Does not actually implement the algorithm (i.e. runBIGEER.R)
  #but the results of those implmentations are included in the inputs directory.

source("scripts//k600_model.R")
rm(list = ls())
print('k model done')

source("scripts//monte_carlo_analysis.R")
rm(list = ls())
print('MC simulations done')

source("scripts//validation_figures.R")
rm(list = ls())
print('Validation figures done. Hopefully you had the correct validation results stored in the output directory!!')
