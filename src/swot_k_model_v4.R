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

#determine if hydraulically wide
data$flag_hydraulicWide <- ifelse(data$Rh/data$depth >= 0.99, 'Rh=H', 'Rh=/=H')

hydraulicallyWide <- filter(data, flag_hydraulicWide == 'Rh=H')

#fit small-eddy model with log-law-of-the-wall dissipation-------------------------------------------------------------------------------------------------------
hydraulicallyWide$logGravitySlope <- log10(g*hydraulicallyWide$slope)
hydraulicallyWide$logDepth <- log10(hydraulicallyWide$depth)
bayeslm_log_hydraulicallyWide_smallEddy_eS <- stan_glm(log10(k600)~logGravitySlope+logDepth, data=hydraulicallyWide,
                                                      family=gaussian(link='identity'),
                                                      prior = normal(location=c(3/8, 1/8), scale=c(1/8,1/8)),
                                                      prior_intercept = normal(location=0, scale=1),
                                                      prior_aux = exponential(rate = 1),
                                                      seed=12345)

#calc 'bayesian r2': https://doi.org/10.1080/00031305.2018.1549100
r2 <- data.frame('r2'=bayes_R2(bayeslm_log_hydraulicallyWide_smallEddy_eS))

#use posterior predicted distributions to predict k600 for plot
predicted_distributions <- as.data.frame(t(posterior_predict(bayeslm_log_hydraulicallyWide_smallEddy_eS)))
predicted_distributions$k600_mean <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = mean, na.rm = T))
predicted_distributions$k600_sigma <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = sd, na.rm = T))
predicted_distributions$k600_05 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.05))}))
predicted_distributions$k600_95 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.95))}))

#add obs k600 to df
predicted_distributions$obs_k600 <- log10(hydraulicallyWide$k600)

#plot
plot_smallEddy_eS <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('Posterior mean\nr2: ', round(mean(r2$r2),2)), x = 0, y = 1.75, size = 8, colour = "#377eb8")+
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = '',
       y = expression(bold(paste(italic(k[600]), '[m/d]'))))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Small-eddy, ', epsilon==epsilon[S]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit reynolds model with log-law-of-the-wall dissipation-------------------------------------------------------------------------------------------------
hydraulicallyWide$logGravitySlope <- log10(g*hydraulicallyWide$slope)
hydraulicallyWide$logDepth <- log10(hydraulicallyWide$depth)
bayeslm_log_hydraulicallyWide_reynolds_eS <- stan_glm(log10(k600)~logGravitySlope+logDepth, data=hydraulicallyWide,
                                                       family=gaussian(link='identity'),
                                                       prior = normal(location=c(9/16, 11/16), scale=c(1/8,1/8)),
                                                       prior_intercept = normal(location=0, scale=1),
                                                       prior_aux = exponential(rate = 1),
                                                       seed=12345)

#calc 'bayesian r2': https://doi.org/10.1080/00031305.2018.1549100
r2 <- data.frame('r2'=bayes_R2(bayeslm_log_hydraulicallyWide_reynolds_eS))

#use posterior predicted distributions to predict k600 for plot
predicted_distributions <- as.data.frame(t(posterior_predict(bayeslm_log_hydraulicallyWide_reynolds_eS)))
predicted_distributions$k600_mean <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = mean, na.rm = T))
predicted_distributions$k600_sigma <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = sd, na.rm = T))
predicted_distributions$k600_05 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.05))}))
predicted_distributions$k600_95 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.95))}))

#add obs k600 to df
predicted_distributions$obs_k600 <- log10(hydraulicallyWide$k600)

#plot
plot_reynolds_eS <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('Posterior mean\nr2: ', round(mean(r2$r2),2)), x = 0, y = 1.75, size = 8, colour = "#377eb8")+
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = '',
       y = '')+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Reynolds extension, ', epsilon==epsilon[S]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit small-eddy model with form-drag dissipation---------------------------------------------------------------------------------
hydraulicallyWide$logGravitySlope <- log10(g*hydraulicallyWide$slope)
hydraulicallyWide$logVms <- log10(hydraulicallyWide$Vms)
bayeslm_log_hydraulicallyWide_smallEddy_eD <- stan_glm(log10(k600)~logGravitySlope+logVms, data=hydraulicallyWide,
                                                      family=gaussian(link='identity'),
                                                      prior = normal(location=c(1/4, 1/4), scale=c(1/8,1/8)),
                                                      prior_intercept = normal(location=0, scale=1),
                                                      prior_aux = exponential(rate = 1),
                                                      seed=12345)

#calc 'bayesian r2': https://doi.org/10.1080/00031305.2018.1549100
r2 <- data.frame('r2'=bayes_R2(bayeslm_log_hydraulicallyWide_smallEddy_eD))

#use posterior predicted distributions to predict k600 for plot
predicted_distributions <- as.data.frame(t(posterior_predict(bayeslm_log_hydraulicallyWide_smallEddy_eD)))
predicted_distributions$k600_mean <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = mean, na.rm = T))
predicted_distributions$k600_sigma <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = sd, na.rm = T))
predicted_distributions$k600_05 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.05))}))
predicted_distributions$k600_95 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.95))}))

#add obs k600 to df
predicted_distributions$obs_k600 <- log10(hydraulicallyWide$k600)

#plot
plot_smallEddy_eD <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('Posterior mean\nr2: ', round(mean(r2$r2),2)), x = 0, y = 1.75, size = 8, colour = "#377eb8")+
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = expression(bold(paste('Predicted ', italic(k[600]), '[m/d]'))),
       y = expression(bold(paste(italic(k[600]), '[m/d]'))))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Small-eddy, ', epsilon==epsilon[D]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit reynolds model with form-drag dissipation------------------------------------------------------------------------------------
hydraulicallyWide$logGravitySlope <- log10(g*hydraulicallyWide$slope)
hydraulicallyWide$logVms <- log10(hydraulicallyWide$Vms)
hydraulicallyWide$logDepth <- log10(hydraulicallyWide$depth)
bayeslm_log_hydraulicallyWide_reynolds_eD <- stan_glm(log10(k600)~logGravitySlope+logVms+logDepth, data=hydraulicallyWide,
                                                      family=gaussian(link='identity'),
                                                      prior = normal(location=c(7/16, 1/4, 9/16), scale=c(1/8, 1/8,1/8)),
                                                      prior_intercept = normal(location=0, scale=1),
                                                      prior_aux = exponential(rate = 1),
                                                      seed=12345)

#calc 'bayesian r2': https://doi.org/10.1080/00031305.2018.1549100
r2 <- data.frame('r2'=bayes_R2(bayeslm_log_hydraulicallyWide_reynolds_eD))

#use posterior predicted distributions to predict k600 for plot
predicted_distributions <- as.data.frame(t(posterior_predict(bayeslm_log_hydraulicallyWide_reynolds_eD)))
predicted_distributions$k600_mean <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = mean, na.rm = T))
predicted_distributions$k600_sigma <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = sd, na.rm = T))
predicted_distributions$k600_05 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.05))}))
predicted_distributions$k600_95 <- (apply(predicted_distributions[1:4000], MARGIN =  1, FUN = function(x){quantile(x, c(0.95))}))

#add obs k600 to df
predicted_distributions$obs_k600 <- log10(hydraulicallyWide$k600)

#plot
plot_reynolds_eD <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('Posterior mean\nr2: ', round(mean(r2$r2),2)), x = 0, y = 1.75, size = 8, colour = "#377eb8")+
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = expression(bold(paste('Predicted ', italic(k[600]), '[m/d]'))),
       y = '')+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Reynolds extension, ', epsilon==epsilon[D]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#combine subplots and write to file
k600_modelPlot_SI <- plot_grid(plot_smallEddy_eS, plot_reynolds_eS, plot_smallEddy_eD, plot_reynolds_eD, ncol=2, label_size = 18, labels=c('a', 'b', 'c', 'd'))
ggsave('cache\\k600_theory\\figS1.jpg', k600_modelPlot_SI, height=9, width=10)

#### SAVE JUST THE FINAL MODEL FOR MAIN TEXT FIGURE---------------------------------------
plot_reynolds_eD <- ggplot(data=predicted_distributions)+
  geom_pointrange(mapping=aes(y=obs_k600, x=k600_mean, xmin=k600_05, xmax=k600_95), fatten=10, color='#377eb8', alpha=0.5) +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = expression(paste('Posterior mean ', r^2, ': 0.50')), x = 0, y = 1.75, size = 8, colour = "#377eb8")+ #hardcoded here to get the right decimals
  scale_y_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  scale_x_continuous(limits=c(-1,2),
                     breaks=c(-1, 0, 1, 2),
                     labels=c('0.1', '1', '10', '100'))+
  labs(x = expression(bold(paste('Predicted ', italic(k[600]), '[m/d]'))),
       y = expression(bold(paste(italic(k[600]), '[m/d]'))))+
  annotation_logticks()+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')
ggsave('cache\\k600_theory\\fig2.jpg', plot_reynolds_eD, height=6, width=6)

#######WRITE MODELS TO FILE-------------------------------
#extract posterior
posterior <- as.data.frame(bayeslm_log_hydraulicallyWide_reynolds_eD)

summary(posterior$logGravitySlope)
summary(posterior$logVms)
summary(posterior$logDepth)
summary(posterior$`(Intercept)`)
summary(posterior$sigma)

models <- data.frame('name'=c('reynolds-eS', 'Small-eddy-eS', 'reynolds-eD', 'Small-eddy-eD'),
                     'r2'=c(mean(bayes_R2(bayeslm_log_hydraulicallyWide_reynolds_eS)),
                            mean(bayes_R2(bayeslm_log_hydraulicallyWide_smallEddy_eS)),
                            mean(bayes_R2(bayeslm_log_hydraulicallyWide_reynolds_eD)),
                            mean(bayes_R2(bayeslm_log_hydraulicallyWide_smallEddy_eD))),
                     'coef'=c(bayeslm_log_hydraulicallyWide_reynolds_eS$coefficients[1],
                              bayeslm_log_hydraulicallyWide_smallEddy_eS$coefficient[1],
                              bayeslm_log_hydraulicallyWide_reynolds_eD$coefficient[1],
                              bayeslm_log_hydraulicallyWide_smallEddy_eD$coefficient[1]),
                     'slope1'=c(bayeslm_log_hydraulicallyWide_reynolds_eS$coefficients[2],
                              bayeslm_log_hydraulicallyWide_smallEddy_eS$coefficient[2],
                              bayeslm_log_hydraulicallyWide_reynolds_eD$coefficient[2],
                              bayeslm_log_hydraulicallyWide_smallEddy_eD$coefficient[2]),
                     'slope2'=c(bayeslm_log_hydraulicallyWide_reynolds_eS$coefficients[3],
                              bayeslm_log_hydraulicallyWide_smallEddy_eS$coefficient[3],
                              bayeslm_log_hydraulicallyWide_reynolds_eD$coefficient[3],
                              bayeslm_log_hydraulicallyWide_smallEddy_eD$coefficient[3]),
                     'slope3'=c(bayeslm_log_hydraulicallyWide_reynolds_eS$coefficients[4],
                              bayeslm_log_hydraulicallyWide_smallEddy_eS$coefficient[4],
                              bayeslm_log_hydraulicallyWide_reynolds_eD$coefficient[4],
                              bayeslm_log_hydraulicallyWide_smallEddy_eD$coefficient[4])
                     )
write.csv(models, 'cache\\k600_theory\\hydraulicWide_models.csv')











#######MONTE CARLO PROPOGATION OF UNCERTANTIES FOR BIKER (i.e. Monte Carol simulation using k600 posterior distributions and Mainngs + rectangular channel uncertainty)-------------------------
posterior <- as.data.frame(bayeslm_log_hydraulicallyWide_reynolds_eD)
u_sigma <- 0.3 #Manning's equation + rectangular channel assumptions

set.seed(13)
n <- 10000 #MC sample size

output <- 1:nrow(hydraulicallyWide)
for (i in output) {
  #bayes regression version
  log_knowns <- sample(posterior$logGravitySlope, n, replace = T)*log(g) + sample(posterior$logGravitySlope,n, replace = T)*log(hydraulicallyWide[i,]$slope) + sample(posterior$logDepth,n, replace = T)*log(hydraulicallyWide[i,]$depth)
  log_k600_pred <- rnorm(n, log(10^(sample(posterior$`(Intercept)`,n, replace = T))) + log_knowns + sample(posterior$logVms,n, replace = T)*rnorm(n, log(hydraulicallyWide[i,]$Vms), u_sigma), sample(log(10^(posterior$sigma)), n, replace=T))
  
  
  output[i] <- sd(log_k600_pred)
}


mean(output)
hist(output)
