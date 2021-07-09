###############
#Set up project
###############
.libPaths('/nas/cee-water/cjgleason/r-lib')
library(ProjectTemplate, lib.loc = "/nas/cee-water/cjgleason/r-lib/")
load.project()

theme_set(theme_cowplot())

#Run analysis-------------------------------
#source(here::here('src' , 'runBIKER.R'))
#source(here::here('src' , 'validation.R'))
source(here::here('src' , 'FCO2_analysis.R'))


#OUT OF DATE I THINK
#source(here::here('src' , 'k600_analysis.R'))
#source(here::here('src' , 'k_eD_comparison.R'))
#source(here::here('src' , 'uncertainty_analysis.R'))
