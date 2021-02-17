##############
#Description: Create maps for manuscript. Can't be run on HPC currently
#Creator: Craig Brinkerhoff
#Date: Winter 2021
##############

#load model settings and packages-------------------------------
source(here::here('scripts' , 'inputs.R'))
library(tmap)
library(st)
library(sf)

#plot measurement locations--------------------------------
tif_sf <- st_as_sf(output, coords = c("lon", "lat"), crs = 4326)
states <- st_read('inputs//MonteCarlo//states.shp')
hydroBox <- st_bbox(tif_sf)

map <- tm_shape(tif_sf, bbox = hydroBox) +
  tm_dots(size=0.1, col = 'darkgreen') +
  tm_scale_bar(position = c('LEFT', 'BOTTOM')) +
  tm_compass(position = c('RIGHT', 'TOP'), size = 1)
rivs <- tm_shape(states, bbox = hydroBox)+
  tm_polygons(col='beige', scale=0.5)
rivs + map
g1 <- grid.grab(width=3, height=3)
