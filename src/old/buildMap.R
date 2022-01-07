#Craig Brinkerhoff
#pepsi 2 river locations map
#Summer 2021

library(ggplot2)
library(tmap)
library(cowplot)
library(RColorBrewer)
library(sf)
library(grid)
theme_set(theme_cowplot())

setwd('C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600')
#setwd("C:/Users/cbrinkerhoff/OneDrive - University of Massachusetts/Ongoing Projects/RSK600")

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
break
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
map <- base + map #+ rivers

tmap_save(map, 'cache/maps/map.jpg')


# 
# 
# 
# #subpanel CA
# ca <- filter(data, River %in% c('AshSlough',
#                                 'Berenda Slough',
#                                 'Chowchilla Canal',
#                                 'Fresno River',
#                                 'Grant Line Canal',
#                                 'Mariposa Bypass',
#                                 'Merced River',
#                                 'Middle River',
#                                 'San Joaquin',
#                                 'San Joaquin 2',
#                                 'Stanislaus River',
#                                 'Tuolumne River',
#                                 'SacramentoDownstream',
#                                 'SacramentoUpstream'))
# box <- st_bbox(ca)
# base <- tm_shape(global_land, bbox = box) +
#   tm_polygons(col = '#99d8c9') +
#   tm_scale_bar(width=1/10, position=c('LEFT', 'BOTTOM'), lwd=2)
# map <- tm_shape(global_rivs, bbox = box)+
#   tm_lines(col='darkblue', scale=0.75)
# rivers <- tm_shape(ca) +
#   tm_symbols(col =  "Data", size=15, shape=21, palette = c('#e31a1c', '#33a02c', '#1f78b4'), border.col = 'black', border.lwd = 1)
# map <- base + map + rivers
# 
# tmap_save(map, 'cache/maps/ca.jpg')
# 
# #---------------------------------------------------------------------------------------------------
# 
# #subpanel eastern US
# es <- filter(data, River %in% c('Iowa River',
#                                 'Mississippi Intermediate',
#                                 'Missouri upstream',
#                                 'Missouri midsection',
#                                 'Missouri downstream',
#                                 'Ohio 1',
#                                 'Ohio 2',
#                                 'Ohio 3',
#                                 'Ohio 4',
#                                 'Ohio 5',
#                                 'Ohio 7',
#                                 'Ohio 8',
#                                 'Olentangy',
#                                 'Connecticut',
#                                 'Cumberland',
#                                 'Kanawha',
#                                 'MississippiDownstream',
#                                 'MississippiUpstream',
#                                 'Ohio',
#                                 'Platte',
#                                 'Wabash'))
# box <- st_bbox(es)
# base <- tm_shape(global_land, bbox = box) +
#   tm_polygons(col = '#99d8c9') +
#   tm_scale_bar(width=1/10, position=c('LEFT', 'BOTTOM'), lwd=2)
# map <- tm_shape(global_rivs, bbox = box)+
#   tm_lines(col='darkblue', scale=0.75)
# rivers <- tm_shape(es) +
#   tm_symbols(col =  "Data", size=15, shape=21, palette = c('#e31a1c', '#33a02c', '#1f78b4'), border.col = 'black', border.lwd = 1)
# map <- base + map + rivers
# 
# tmap_save(map, 'cache/maps/es.jpg')
# 
# #-------------------------------------------------------------------------------
# 
# #subpanel europe
# eu <- filter(data, River %in% c('Seine Downstream',
#                                 'Seine Upstream',
#                                 'GaronneDownstream',
#                                 'Po',
#                                 'Seine',
#                                 'Severn'))
# box <- st_bbox(eu)
# box[1] <- box[1] - 1
# base <- tm_shape(global_land, bbox = box) +
#   tm_polygons(col = '#99d8c9') +
#   tm_scale_bar(width=1/10, position=c('LEFT', 'BOTTOM'), lwd=2)
# map <- tm_shape(global_rivs, bbox = box)+
#   tm_lines(col='darkblue', scale=0.75)
# rivers <- tm_shape(eu) +
#   tm_symbols(col =  "Data", size=15, shape=21, palette = c('#e31a1c', '#33a02c', '#1f78b4'), border.col = 'black', border.lwd = 1)
# map <- base + map + rivers
# 
# tmap_save(map, 'cache/maps/eu.jpg')
# 
# #-------------------------------------------------------------------------------
# 
# #subpanel bangladesh
# bg <- filter(data, River %in% c('Brahmaputra',
#                                 'Jamuna',
#                                 'Kushiyara',
#                                 'Padma',
#                                 'Ganges'))
# box <- st_bbox(bg)
# base <- tm_shape(global_land, bbox = box) +
#   tm_polygons(col = '#99d8c9') +
#   tm_scale_bar(width=1/10, position=c('LEFT', 'BOTTOM'), lwd=2)
# map <- tm_shape(global_rivs, bbox = box)+
#   tm_lines(col='darkblue', scale=0.75)
# rivers <- tm_shape(bg) +
#   tm_symbols(col =  "Data", size=15, shape=21, palette = c('#e31a1c', '#33a02c', '#1f78b4'), border.col = 'black', border.lwd = 1)
# map <- base + map + rivers
# 
# tmap_save(map, 'cache/maps/bg.jpg')
