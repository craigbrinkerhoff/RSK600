##################
## Calculates SWORD/SWOT mean annual surface area (widths come from GRWL and /or MERIT Hydro) and network length
## Does each continent separately to spare holding the global river network in memory...
## Craig Brinkerhoff
## Summer 2022
###################

library(sf)

#path to SWORD v6
setwd('C:\\Users\\craig\\OneDrive - University of Massachusetts\\Datasets\\SWORD_shp_v2\\shp')

########################################SWOT/SWORD########################################################################
#get files
files <- list.files(recursive = TRUE)
files <- files[substr(files, nchar(files)-3, nchar(files)) == '.shp']
files <- files[substr(files, 13, 19) == 'reaches']

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
  
  sword <- sword[sword$width >= 50 & substr(sword$reach_id, 11,11) == '1',] #grab only rivers, ignore
  
  #estimate SA
  SA <- sum((sword$width * sword$reach_len)*1e-6, na.rm=T) #[m2]
  #SA <- SA/1000000 #[km2]
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
length_fin/443509286 #Liu et al 2022

df <- data.frame('sa_km2'=c(SA_fin, 811000),
                 'length_km' = c(length_fin, 443509286),
                 'model'=c('swot', 'total'))

write.csv(df, 'C:/Users/craig/OneDrive - University of Massachusetts/Ongoing Projects/RSK600/cache/swot_total_comparison.csv')
