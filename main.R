#Runs the three analysis scripts for this project.
#Does not actually implement the algorithm (i.e. runBIGEER.R)
  #but the results of those implmentations are included in the inputs directory.

#k600 model development
source("scripts//k600_model.R")
rm(list = ls())
print('k model done')

#Uncertainity analysis
source("scripts//monte_carlo_analysis.R")
rm(list = ls())
print('MC simulations done')

#Produce validation figures
source("scripts//validation_figures.R")
rm(list = ls())
print('Validation figures done. Hopefully you had the correct validation results stored in the output directory!!')

#Move all results to global environment
results_3_1 <- read.csv('outputs//k600//results.csv')
results_3_1_MC <- read.csv('outputs//MonteCarlo//results.csv')
results_3_2_rivs <- read.csv('outputs//validation//results_by_riv.csv')
results_3_2_all <- read.csv('outputs//validation//results_all_riv.csv')
results_3_3_all <- read.csv('outputs//flux_implications//fco2_stats_all.csv')
results_3_3_bulk <- read.csv('outputs//flux_implications//bulkFluxes.csv')
results_3_3_bulk$key <- as.character(results_3_3_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v1.Rmd")
print('Project report, or Outline.html, is generated.')