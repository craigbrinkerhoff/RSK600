#################
#This will grab and store all results in your current RStudio environment so that you can knit the manuscript correctly. It then knits the manuscript to a word doc
#################

#Move all results to global environment
modelsWIDE <- read.csv('cache\\k600_theory\\hydraulicWide_models.csv')
n <- 166 #number of hydraulically wide measurements used in k600 validation
HG_swot <- read.csv('cache\\k600_theory\\HG_swot.csv')

results_k600_rivs <- read.csv('cache//validation//results_by_riv.csv')
results_fco2_bulk <- read.csv('cache//FCO2//massFluxes.csv')
results_fco2_bulk$key <- as.character(results_fco2_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_r_v1.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_r_v1.Rmd")
print('Supplemental Information generated.')
