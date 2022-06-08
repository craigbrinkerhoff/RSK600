#################
##Description: Development of gas exchange velocity model for 'hydraulically wide channels'
              #It also estimates how many SWOT-observable rivers meet this hydraulically-wide condition
##Creator: Craig Brinkerhoff
##Date: Fall 2021

##TO RUN:
#  1) set the working directory below to the '~/RSK600' directory
#  2) run script!
#################

#######SETUP-----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(colorspace)
library(cowplot)
library(readr)
library(rstanarm)
library(grid)
theme_set(theme_classic())

#some constants
g <- 9.8 #gravitational acceleration [m/s2]

#######READ IN HYDRAULICS DATA FROM BRINKERHOFF ETAL 2019 TO GET LARGER HYDRAULICS DATASET------------------------------
data <- read.csv('data\\Brinkerhoff_etal_2019\\field_measurements.csv')

#some necessary filtering
data <- filter(data, is.finite(chan_width) ==1) %>%
  filter(is.finite(chan_velocity)==1) %>%
  filter(is.finite(chan_discharge)==1) %>%
  filter(is.finite(chan_area)==1) %>%
  filter(chan_width > 0) %>%
  filter(chan_velocity > 0) %>%
  filter(chan_discharge > 0) %>%
  filter(chan_area > 0) %>%
  filter(measured_rating_diff %in% c('Excellent', 'EXCL', 'GOOD', 'Good'))

#imperial to metric
data$area <- data$chan_area * 0.092903 #ft2 to m2
data$width <- data$chan_width*0.305 #m
data$Vms <- data$chan_velocity*0.305 #m/s
data$depth <- data$area / data$width
data$slope <- NA 
data$Qm3s <- data$chan_discharge * 0.0283 #ft3/s to m3/s

#calculate bankfull width so we can throw out out-of-bank events
bankfullW = group_by(data, site_no) %>%
  filter(n() >= 20) %>%
  mutate(rank = rank(width, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_width = max(width[temp]))
data <- left_join(data, bankfullW, 'site_no')

#calculate bankfull depth so we can throw out out-of-bank events
bankfullD = group_by(data, site_no) %>%
  filter(n() >= 20) %>%
  mutate(rank = rank(depth, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_depth = max(depth[temp]))
data <- left_join(data, bankfullD, 'site_no')

data <- filter(data, width < bankful_width)
data <- filter(data, depth < bankful_depth)

#######READ IN FIELD DATA COLLECTED FROM USGS REPORT-------------------
  #Chuchill et al 1964
  #Owens et al 1964
usgs_data <- read.csv('data/additional_USGS_reaeration_data.csv')
usgs_data$depth <- usgs_data$depth_ft * 0.3048 #ft to m
usgs_data$width <- usgs_data$width_ft * 0.3048 #ft to m
usgs_data$Vms <- usgs_data$velocity_fps * 0.3048 #ft to m
usgs_data$slope <- usgs_data$slope_10_x_4 * 10^-4
usgs_data$Qm3s <- usgs_data$depth * usgs_data$width * usgs_data$Vms #[m3/s]
usgs_data$k600 <- usgs_data$k20_1_day * usgs_data$depth * (600/530)^(-1/2) #convert to k in m/dy and to a Sc of 600
usgs_data_s <- select(usgs_data, 'width', 'Vms', 'depth', 'slope', 'k600', 'Qm3s', 'study')

#######READ IN FIELD DATA WITH K600-------------------
  #Ulseth etal 2019
ulseth_data <- read.csv('data/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM")
ulseth_data <- ulseth_data[,-8]
ulseth_data$dataset <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$dataset == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$dataset)) #contains data from 5 studies
ulseth_data <- filter(ulseth_data, is.na(width)==0) #filter out no width measurements
ulseth_data_s <- select(ulseth_data, 'width', 'Vms', 'depth', 'slope', 'k600', 'Qm3s')
ulseth_data_s$study <- 'Ulseth_etal_2019'

#######JOIN DATASETS---------------------------------
data <- select(data, 'width', 'Vms', 'depth', 'slope', 'Qm3s')
data$k600 <- NA
data$study <- 'Brinkerhoff_etal_2019'
data <- rbind(data, ulseth_data_s, usgs_data_s)

#########QUICK DETOUR TO CALCULATE HG MODELS FOR FCO2 CALCULATIONS (DURING BIKER VALIDATION)------------------------------
num <- nrow(data)
wid <- lm(log(width)~log(Qm3s), data=data)
dep <- lm(log(depth)~log(Qm3s), data=data)
write_rds(wid, 'cache/widAHG.rds')
write_rds(dep, 'cache/depAHG.rds')

#######CALCULATE CHANNEL GEOMETRY---------------------------------------------------------
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius, assuming rectangular channel [m]
data$ustar <- sqrt(g*data$Rh*data$slope) #friction/shear velocity [m/s]

#######SOME FLAGS FOR SWOT-OBSERVABLE RIVERS AND WHEN RH=H-------------------------------
data$flag_swot <- ifelse(data$width >= 100, 'SWOT', 'Small')

#########HYDRAULIC GEOMETRY OF SWOT-OBSERVABLE FLOWS
HG_swot_mean <- round(mean((data[data$flag_swot == 'SWOT',]$Rh/data[data$flag_swot == 'SWOT',]$depth)),2)
HG_swot_sd <- round(sd((data[data$flag_swot == 'SWOT',]$Rh/data[data$flag_swot == 'SWOT',]$depth)),2)
n <- nrow(data[data$flag_swot == 'SWOT',])
df <- data.frame('mean'=HG_swot_mean, 'sd'=HG_swot_sd, 'n'=n)
write.csv(df, 'cache/k600_theory/HG_swot.csv')

#######LOG TRANSFORM SOME VARIABLES-----------------------------------------------------------
data$log_slope <- log(data$slope)

######## K600 MODELS--------------------
#hydraulically-wide rivers
data <- filter(data, study %in% c('Ulseth_etal_2019', 'Churchill_etal_1962', 'Owens_etal_1964'))
out <- data.frame('width_min'=min(data$width),
                  'width_max'=max(data$width),
                  'width_mean'=mean(data$width),
                  'Q_min'=min(data$Qm3s),
                  'Q_max'=max(data$Qm3s),
                  'Q_mean'=mean(data$Qm3s),
                  'k600_min'=min(data$k600),
                  'k600_max'=max(data$k600),
                  'k600_mean'=mean(data$k600))
write.csv(out, 'cache/k600_dataStats.csv') #save dataset stats to file

#determine if hydraulically wide
data$flag_hydraulicWide <- ifelse(data$Rh/data$depth >= 0.99, 'Rh=H', 'Rh=/=H')

hydraulicallyWide <- filter(data, flag_hydraulicWide == 'Rh=H')

#tests for reviewer 1############################################
#test 1
hydraulicallyWide$log_hydraulicallyWideModel <- (7/16)*log(g*hydraulicallyWide$slope) + (1/4)*log(hydraulicallyWide$Vms) + (9/16)*log(hydraulicallyWide$depth)
beta_1 <- mean(c(log(hydraulicallyWide$k600) - hydraulicallyWide$log_hydraulicallyWideModel))
hydraulicallyWide$logk600_pred_wideHydraulics <- beta_1+hydraulicallyWide$log_hydraulicallyWideModel
lm_log_hydraulicallyWide_reynolds_eD <- lm(log(k600)~logk600_pred_wideHydraulics, data=hydraulicallyWide)

summary(lm_log_hydraulicallyWide_reynolds_eD)

test1Plot <- ggplot(hydraulicallyWide, aes(y=k600, x=exp(logk600_pred_wideHydraulics)))+
  geom_point(size=5, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_log_hydraulicallyWide_reynolds_eD)$r.squared,2)), x = 1, y = 100, size = 8, colour = "#377eb8")+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  labs(x = '',
       y = 'Observed k600')+
  annotation_logticks()+
  ggtitle(expression(bold('Test 1')))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#test 2
hydraulicallyWide$term1 <- log(g*hydraulicallyWide$slope)
hydraulicallyWide$term2 <- log(hydraulicallyWide$Vms)
hydraulicallyWide$term3 <- log(hydraulicallyWide$depth)
lm_log_hydraulicallyWide_reynolds_eD <- lm(log(k600)~term1+term2+term3, data=hydraulicallyWide)
hydraulicallyWide$logk600_pred_wideHydraulics <- predict(lm_log_hydraulicallyWide_reynolds_eD, hydraulicallyWide)

summary(lm_log_hydraulicallyWide_reynolds_eD)

test2Plot <- ggplot(hydraulicallyWide, aes(y=k600, x=exp(logk600_pred_wideHydraulics)))+
  geom_point(size=5, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_log_hydraulicallyWide_reynolds_eD)$r.squared,2)), x = 1, y = 100, size = 8, colour = "#377eb8")+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  labs(x = 'Predicted k600',
       y = '')+
  annotation_logticks()+
  ggtitle(expression(bold('Test 2')))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#test 3
hydraulicallyWide$logGravitySlope <- log10(g*hydraulicallyWide$slope)
hydraulicallyWide$logVms <- log10(hydraulicallyWide$Vms)
hydraulicallyWide$logDepth <- log10(hydraulicallyWide$depth)
bayeslm_log_hydraulicallyWide_reynolds_eD <- stan_glm(log10(k600)~logGravitySlope+logVms+logDepth, data=hydraulicallyWide,
                                                      family=gaussian(link='identity'),
                                                      prior = normal(location=c(7/16, 1/4, 9/16), scale=c(1/8, 1/8,1/8)),
                                                      prior_intercept = normal(location=0.47, scale=1),#location=0 after internal recentering
                                                      #intercept prior is left weakly informed using default prior, i.e. normal(0.47, 2.5)
                                                      #sigma prior is left weakly informed using default prior, i.e. exponential(rate=1)
                                                      seed=12345)

#calc 'bayesian r2': https://doi.org/10.1080/00031305.2018.1549100
r2 <- data.frame('r2'=bayes_R2(bayeslm_log_hydraulicallyWide_reynolds_eD))
r2Plot <- ggplot(r2, aes(x=r2))+
  geom_density(aes(y=..scaled..),fill='skyblue4', color='black', size=1.25) +
  geom_vline(xintercept = median(r2$r2), linetype='dotted', size=2)+
  ylab('Scaled Density')+
  xlab('Posterior r2')+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#posterior predictive check to give a 'Bayesian' check
posteriorCheckPlot <- pp_check(bayeslm_log_hydraulicallyWide_reynolds_eD, nreps=100)+
  xlab('log10 k600 [m/dy]') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = c(.9, .9))

postVpriorPlot <- posterior_vs_prior(bayeslm_log_hydraulicallyWide_reynolds_eD) +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#extract posterior
posterior <- as.data.frame(bayeslm_log_hydraulicallyWide_reynolds_eD)

summary(posterior$logGravitySlope)
summary(posterior$logVms)
summary(posterior$logDepth)
summary(posterior$`(Intercept)`)
summary(posterior$sigma)

#use posterior predicted distributions to predict k600 for plot
predicted_distributions <- as.data.frame(t(posterior_predict(bayeslm_log_hydraulicallyWide_reynolds_eD)))
predicted_distributions$k600_mean <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = mean, na.rm = T))
predicted_distributions$k600_sigma <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = sd, na.rm = T))
predicted_distributions$k600_05 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.05))}))
predicted_distributions$k600_95 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.95))}))

#add obs k600 to df
predicted_distributions$obs_k600 <- log10(hydraulicallyWide$k600)

#plot
test3Plot <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('Posterior mean\nr2: ', round(mean(r2$r2),2)), x = 0, y = 1.75, size = 8, colour = "#377eb8")+
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = 'Predicted k600',
       y = 'Observed k600')+
  annotation_logticks()+
  ggtitle(expression(bold('Test 3')))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#text description
description <- textGrob('Test 1: linear regression\nusing the logged equation 7\n\nTest 2: linear regression\nusing the logged equation 7 terms,\nbut the exponents can vary\n\nTest 3: Bayesian linear regression\nusing the logged equation 7,\nwhere the theoretically-defensible\nexponents are set as the priors', y=0.5, gp=gpar(fontface="bold", col="darkblue", fontsize=18))


#bring it all together
grid <- plot_grid(test1Plot, test2Plot, test3Plot, description, ncol=2)
ggsave('docs/review_round2/modelComparsion.jpg', grid, width=10, height=10)

bayesGrid <- plot_grid(r2Plot, posteriorCheckPlot, ncol=2)
bayesGrid <- plot_grid(bayesGrid, postVpriorPlot, ncol=1)
ggsave('docs/review_round2/bayesModelCheck.jpg', bayesGrid, width=10, height=10)
