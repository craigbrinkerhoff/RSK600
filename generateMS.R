#################
#This will generate and store all results in your current RStudio environment so that you can knit the manuscript correctly, and then it'' knit the manuscript too.
#################

setwd('C:\\Users\\craig\\Documents\\GitHub\\RSK600')

#Move all results to global environment
#results_3_1 <- read.csv('cache//k600//results.csv')
#results_3_1_MC <- read.csv('cache//MonteCarlo//results.csv')
results_3_2_rivs <- read.csv('cache//validation//results_by_riv.csv')
results_3_2_all <- read.csv('cache//validation//results_all_riv.csv')
results_3_3_all <- read.csv('cache//FCO2//fco2_stats_all.csv')
results_3_3_rivs <- read.csv('cache//FCO2//fco2_stats_by_river.csv')
results_3_3_bulk <- read.csv('cache//FCO2//bulkFluxes.csv')
results_3_3_bulk$key <- as.character(results_3_3_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v3.Rmd")
print('Manuscript generated.')

rmarkdown::render("docs//supp_v1.Rmd")
print('Supplemental Information generated.')
