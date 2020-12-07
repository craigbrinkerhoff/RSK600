#Runs the three analysis scripts for this project.
#Does not actually implement the algorithm (i.e. runBIGEER.R)
  #but the results of those implmentations are included in the inputs directory.

source("scripts//k600_model.R")
rm(list = ls())
print('k model done')

source("scripts//monte_carlo_analysis.R")
rm(list = ls())
print('MC simulations done')

source("scripts//Figures_pepsi.R")
rm(list = ls())
print('Validation figures done. Hopefully you had the validation results stored in the project!!')
