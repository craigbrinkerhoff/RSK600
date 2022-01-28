#################
#This will grab and store all results in your current RStudio environment so that you can knit the manuscript correctly. It then knits the manuscript to a
  #word doc that, via the redoc package, can be tracked changed and convert back to Rmd!
#################

#Move all results to global environment
modelsWIDE <- read.csv('cache\\k600_theory\\hydraulicWide_models.csv')
n <- 166 #number of hydraully wide measurements used in k600 validation
HG_swot <- read.csv('cache\\k600_theory\\HG_swot.csv')

results_dynamics <- read.csv('cache/validation/results_dynamics.csv')
results_k600_rivs <- read.csv('cache//validation//results_by_riv.csv')
results_fco2_rivs <- read.csv('cache//FCO2//results_by_riv.csv')
results_fco2_bulk <- read.csv('cache//FCO2//massFluxes.csv')
results_fco2_bulk$key <- as.character(results_fco2_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v8.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_v8.Rmd")
print('Supplemental Information generated.')
