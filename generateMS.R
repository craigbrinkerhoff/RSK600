#################
#This will grab and store all results in your current RStudio environment so that you can knit the manuscript correctly. It then knits the manuscript too.
#################

library(dplyr)

#setwd("C:/Users/cbrinkerhoff/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")
setwd("C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")

#Move all results to global environment
modelsWIDE <- read.csv('cache\\k600_theory\\hydraulicWide_models.csv')
percs <- read.csv('cache\\k600_theory\\Rh_H_percents.csv')
percs_overall <- group_by(percs, flag_hydraulicWide) %>% 
  summarise(total = sum(n)) %>%
  mutate(perc = 100 * total/sum(total))

results_k600_rivs <- read.csv('cache//validation//results_by_riv.csv')
results_fco2_rivs <- read.csv('cache//FCO2//CO2_results.csv')
results_fco2_bulk <- read.csv('cache//FCO2//massFluxes.csv')
results_fco2_bulk$key <- as.character(results_fco2_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v6.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_v6.Rmd")
print('Supplemental Information generated.')


#scraps
#results_k600_all <- read.csv('cache//validation//results_all_riv.csv')
#results_fco2_all <- read.csv('cache//FCO2//fco2_stats_all.csv')