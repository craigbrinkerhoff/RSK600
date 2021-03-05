#################
#uncomment this code if you want to run the project on your local machine. Optimized for parallel computing on an HPC so only the strongest machines will survive!!!
#################


# #k600 model development
# source("scripts//k600_analysis.R")
# rm(list = ls())
# 
# #Uncertainity analysis
# source("scripts//MonteCarlo_analysis.R")
# rm(list = ls())
# 
# #run BIKER
# source("scripts//runBIKER.R")
# rm(list = ls())
# 
# #Produce validation figures
# source("scripts//validation.R")
# rm(list = ls())
# 
# #FCo2 analysis
# source("scripts//FCO2_analysis.R")
# rm(list = ls())

#Move all results to global environment
results_3_1 <- read.csv('outputs//k600//results.csv')
results_3_1_MC <- read.csv('outputs//MonteCarlo//results.csv')
results_3_2_rivs <- read.csv('outputs//validation//results_by_riv.csv')
results_3_2_all <- read.csv('outputs//validation//results_all_riv.csv')
results_3_3_all <- read.csv('outputs//FCO2//fco2_stats_all.csv')
results_3_3_rivs <- read.csv('outputs//FCO2//fco2_stats_by_river.csv')
results_3_3_bulk <- read.csv('outputs//FCO2//bulkFluxes.csv')
results_3_3_bulk$key <- as.character(results_3_3_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v2.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_v1.Rmd")
print('Supplemental Information generated.')