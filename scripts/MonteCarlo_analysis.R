#Monte Carlo analysis for quantifying uncertainity of CO2 evasion model
  #This quantifies model error for two scenarios: 1) uncertainity associated with Manning's equation, and
    #2) uncertainity associated with Manning's equation AND the gas transfer velocity model from Raymond et al (2012)
    #The former is used when we are assuming that there is no error in the gas transfer velocity equation.
    #We do this to calculate Flux using the same model to test the use of a Bayesian framework
#Creator: Craig Brinkerhoff
#Summer 2020

library(cowplot)
library(tidyverse)
library(tmap)
library(st)
library(sf)
library(grid)
theme_set(theme_cowplot())

#random samplers: M sets of reasonable observations for model
M <- 8000 # number of sets of hydraulic 'measurements'
set.seed(455)

#Pull in hydraulic measurements from Brinkerhoff etal 2019 to test parameter uncertainity------------------------------------------
hydraulics_init <- read.csv('/home/cbrinkerhoff_umass/edu/datasets/MonteCarlo/field_measurements.csv')
nhd_med <- read.csv('inputs//MonteCarlo//NHD_join_table.csv')

#add reach slopes from NHD
hydraulics <- left_join(hydraulics_init, select(nhd_med, SLOPE, SOURCE_FEA, LatSite, LonSite), by = c("site_no" = "SOURCE_FEA"))

#intial cleaning of field measurements
hydraulics <- filter(hydraulics, is.finite(chan_width) ==1) %>%
  filter(is.finite(chan_velocity)==1) %>%
  filter(is.finite(chan_discharge)==1) %>%
  filter(is.finite(chan_area)==1) %>%
  filter(chan_width > 0) %>%
  filter(chan_velocity > 0) %>%
  filter(chan_discharge > 0) %>%
  filter(chan_area > 0) %>%
  filter(measured_rating_diff != 'Poor') %>%
  filter(measured_rating_diff != 'POOR') %>%
  filter(is.finite(SLOPE)==1)

#imperial to metric
hydraulics$area <- hydraulics$chan_area * 0.092903 #ft2 to m2
hydraulics$width <- hydraulics$chan_width*0.305 #m
hydraulics$velocity <- hydraulics$chan_velocity*0.305 #m/s
hydraulics$slope <- hydraulics$SLOPE
hydraulics$depth <- hydraulics$area / hydraulics$width #m
hydraulics$n <- (hydraulics$depth^(2/3)*hydraulics$SLOPE^(1/2))/hydraulics$velocity

#many slopes were saved with floor values of 1e-5, so I just remove these because they influence the k600 model to provide artifical results
hydraulics <- filter(hydraulics, slope > 1e-5)

#randomly grab M hydraulic sets from the Brinkerhoff etal 2019 dataset (which is 500,000+...)
ids <- runif(M, 1, nrow(hydraulics))
slope <- hydraulics$slope[ids]
velocity <- hydraulics$velocity[ids]
width <- hydraulics$width[ids]
lats <- hydraulics$LatSite[ids]
lons <- hydraulics$LonSite[ids]

#Model Parameters-----------------------------
man_sd <- 0.25 #Hagemann et al. 2017's manning's equation uncertainty ("from HydroSWOT")

#Monte Carlo simulations (M simulations of 10,000 runs each)--------------------------------------------------------------
output <- data.frame('sigma_logk600'=NA, 'expectation_logk600'=NA, 'vel'=NA, 'slope'=NA, 'width'=NA, 'lat'=NA, 'lon'=NA)
for(i in 1:M){
  vel_man <- exp(rnorm(10000, log(velocity[i]), man_sd))
  
  if(width < 10 & slope < 0.05){ #< 10m wide and slope < 0.05
    exp <- rnorm(10000, 0.6488131, 0.173) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 111.58121, 1.038) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  if(width < 10 & slope >= 0.05){ #< 10m wide and slope >= 0.05
    exp <- rnorm(10000, 1.3160065, 0.177) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 792.63149, 1.087) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  if(width < 50 & width >= 10){ #10-50m wide
    exp <- rnorm(10000, 0.6624354, 0.241) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 109.04977, 1.046) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  if(width < 100 & width >= 50){ #50-100m wide
    exp <- rnorm(10000, 0.4366114, 0.359) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 31.84344, 1.059) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  if(width >= 100){ #> 100m wide
    exp <- rnorm(10000, 0.2724834, 0.201) #Sample from normal distribution. sigma taken from my rule-based regression model
    int <- rnorm(10000, 14.16939, 1.031) #Sample from normal distribution. sigma taken from my rule-based regression model
  }
  
  #k600 model
  log_eD <- log(slope[i] * 9.8 * vel_man)
  logk <- exp*log_eD + log(int)
  
  #get model estimate uncertainty
  temp <- c(sd(logk, na.rm=T), mean(logk, na.rm = T), velocity[i], slope[i], width[i], lats[i], lons[i])
  output <- rbind(output, temp)
  
  #save a few random uncertantity scenarios
  if (i == 100) {
    logk_100 <- logk
  }
  if (i == 3700) {
    logk_3700 <- logk
  }
  if (i == 4000) {
    logk_4000 <- logk
  }
}

output <- output[-1,]

#Plot histogram of M uncertanties--------------------------
meanSD <- mean(output$sigma_logk600, na.rm=T)
medianSD <- median(output$sigma_logk600, na.rm=T)
plot <- ggplot(output, aes(x=sigma_logk600)) +
  geom_histogram(color='black', fill='darkgreen', size=1, bins=50) +
  xlab('lnSD of MC Simulated Estimates [m/dy]') +
  ylab('Count') +
  geom_vline(xintercept = meanSD, linetype='dashed', size=1.2, color='blue')

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

#Plot some example uncertanties--------------------------------
t <- data.frame(logk_100, logk_3700, logk_4000)
t <- gather(t, key=key, value=value)
plot2 <- ggplot(t, aes(x=value, color=key)) +
  geom_density(size=2) +
  ylab('Density') +
  scale_color_brewer(palette = 'Dark2') +
  xlab('ln k600 [m/dy]') +
  theme(legend.position = 'none')

combinedPlot <- plot_grid(plot, plot2, ncol=2, labels = 'auto', label_size = 18)
combinedPlot <- plot_grid(combinedPlot, g1, ncol=1, labels = c(NA, 'c'), label_size = 18)

#combine subplots-------------------------------------------------
ggsave('outputs//MonteCarlo//MCsimulations.jpg', combinedPlot, width=10, height=6)

#Save uncertainity to file for ms---------------------------------------
results <- data.frame('meanSigma'=meanSD, 'medianSigma'=medianSD)
write.csv(results, 'outputs//MonteCarlo//results.csv')
