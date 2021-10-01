###############
#Set up project
###############
.libPaths('/nas/cee-water/cjgleason/r-lib')
library(ProjectTemplate, lib.loc = "/nas/cee-water/cjgleason/r-lib/")
load.project()

theme_set(theme_cowplot())

#Run analysis-------------------------------
source(here::here('src' , 'runBIKER.R'))
source(here::here('src' , 'FCO2_analysis.R'))
source(here::here('src' , 'figures_k600.R'))
source(here::here('src' , 'figures_FCO2.R'))
