##################
## Calculates SWORD/SWOT mean annual surface area (widths come from GRWL and /or MERIT Hydro) and network length
## Does eachcontinent seperately to spare holding the global river network in memory...
## Craig Brinkerhoff
## Summer 2022
###################

library(sf)

#path to SWORD v6
setwd('C:\\Users\\cbrinkerhoff\\OneDrive - University of Massachusetts\\Datasets\\SWORD_v06\\Reaches_Nodes\\shp')

########################################SWOT/SWORD########################################################################
#get files
files <- list.files(pattern = '*.shp', recursive = TRUE)
files <- files[substr(files, 22, 28) == 'reaches']

#run per continent, concatenating SA
SA_fin <- 0
length_fin <- 0
continent <- unique(substr(files, 1, 2))
for (k in continent){
  files2 <- files[substr(files, 1, 2) == k]
  
  #concatenate shapefiles
  sword <- st_read(files2[1], quiet=TRUE)
  files2 <- files2[-1]
  for(i in files2){
    temp_sf <- st_read(i, quiet=TRUE)
    sword <- rbind(sword, temp_sf)
  }
  
  #estimate SA
  SA <- sum(sword$width * sword$reach_len, na.rm=T) #[m2]
  SA <- SA*1e-6 #[km2]
  SA_fin <- SA_fin + SA
  print(paste0(k, ' SA: ', SA))
  print(paste0('total SA: ', SA_fin))
  
  #estimate length
  length <- sum(sword$reach_len) #[m]
  length <- length / 1000 #[km]
  length_fin <- length + length_fin
  print(paste0(k, ' Length: ', length))
  print(paste0('total Length: ', length_fin))
}

#calc percent SWOT observable
SA_fin/811000 #Liu et al 2022
length_fin/23300000 #riverATLAS Messager et al 2021
