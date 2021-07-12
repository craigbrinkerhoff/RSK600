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
