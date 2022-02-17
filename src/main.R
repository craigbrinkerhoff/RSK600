###############
#Set up and run project
###############
.libPaths('/nas/cee-water/cjgleason/r-lib')
library(ProjectTemplate, lib.loc = "/nas/cee-water/cjgleason/r-lib/")
load.project()

theme_set(theme_cowplot())
start <- Sys.time()

#Run analysis-------------------------------
source(here::here('src' , 'runBIKER.R'))
source(here::here('src' , 'runFCO2.R'))
source(here::here('src' , 'figures_k600.R'))
source(here::here('src' , 'figures_FCO2.R'))

end <- Sys.time()
end - start
