#################
#This will grab and store all results in your current RStudio environment so that you can knit the manuscript correctly. It then knits the manuscript too.
#################

setwd("C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")

#Move all results to global environment
#turbulence_stats <- read.csv('cache\\k600_theory\\turbulence_stats.csv')
models_ulseth <- read.csv('cache\\k600_theory\\ulseth_models.csv')
models_ustar <- read.csv('cache\\k600_theory\\ustar_models.csv')
models_eD4 <- read.csv('cache\\k600_theory\\eD_models.csv')
percs <- read.csv('cache\\k600_theory\\Rh_H_percents.csv')
percs_overall <- dplyr::group_by(percs, flag_depth) %>% 
  dplyr::summarise(total = sum(n)) %>%
  dplyr::mutate(perc = 100 * total/sum(total))

results_3_2_rivs <- read.csv('cache//validation//results_by_riv.csv')
results_3_2_all <- read.csv('cache//validation//results_all_riv.csv')
results_3_3_all <- read.csv('cache//FCO2//fco2_stats_all.csv')
results_3_3_bulk <- read.csv('cache//FCO2//bulkFluxes.csv')
results_3_3_bulk$key <- as.character(results_3_3_bulk$key)

#knit manuscript
rmarkdown::render("manuscript_v5.Rmd")
print('Manuscript generated.')

rmarkdown::render("supp_v5.Rmd")
print('Supplemental Information generated.')
