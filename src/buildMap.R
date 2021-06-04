#Craig Brinkerhoff
#pepsi 2 river locations map
#Summer 2021

library(tidyverse)
library(tmap)
library(cowplot)
library(RColorBrewer)
library(sf)
library(grid)
theme_set(theme_cowplot())

setwd('C:\\Users\\craig\\Documents\\GitHub\\RSK600')

##########
#Read in locations and rivers----------------
#########
data <- read.csv('data\\dataLocations.csv', header = T)
colnames(data) <- c('River', 'lat', 'lon', 'Data')
#data <- data[-1,]
data$lat <- as.numeric(as.character(data$lat))
data$lon <- as.numeric(as.character(data$lon))

global_rivs <- st_read('data\\ne_10m_rivers_lake_centerlines\\ne_10m_rivers_lake_centerlines.shp')
global_land <- st_read('data\\ne_50m_land\\ne_50m_land.shp')
data <- st_as_sf(data, coords = c('lon', 'lat'))
st_crs(data) <- st_crs(global_rivs)

##########
#main map------------
##########
base <- tm_shape(global_land) +
  tm_polygons(col = '#99d8c9') +
  tm_compass(position = c('RIGHT', 'TOP')) +
  tm_scale_bar(width=1/10, position=c('RIGHT', 'BOTTOM'), lwd=2)
map <- tm_shape(global_rivs)+
  tm_lines(col='darkblue', scale=0.75)
rivers <- tm_shape(data) +
  tm_symbols(col =  "Data", size=0.5, shape=21, palette = c('#e31a1c', '#33a02c', '#1f78b4'), border.col = 'black', border.lwd = 1)
map <- base + map + rivers

tmap_save(map, 'cache/map.jpg')