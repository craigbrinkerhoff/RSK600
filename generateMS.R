#################
#This will grab and store all results in your current RStudio environment so that you can knit the manuscript correctly. It then knits the manuscript to a
  #word doc that, via the redoc package, can be tracked changed and convert back to Rmd!
#################

#setwd("C:/Users/cbrinkerhoff/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")
setwd("C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")

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
rmarkdown::render("manuscript_v7.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_v7.Rmd")
print('Supplemental Information generated.')


#To convert tracked changes edited manuscripts back to rmarkdown! #https://noamross.github.io/redoc/articles/mixed-workflows-with-redoc.html
#NOTES:
  #Can also be set to automatically accept tracked changes and comments.
  #This package is no longer maintained...
  #for tracked changes:
    #{-- --}
  #comments:
    #comments in {== highlighted text  ==}{>> comment <<}

#redoc::dedoc('manuscript_v7_e.docx', track_changes = 'criticmarkup', overwrite=TRUE, wrap=getOption("redoc.wrap", NA))
