#################
##Description: Developing gas exchange model for 'hydraulically wide channels'
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
theme_set(theme_classic())

#setwd('C:\\Users\\cbrinkerhoff\\OneDrive - University of Massachusetts\\Ongoing Projects\\RSK600')
setwd('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Ongoing Projects\\RSK600')

#some constants
g <- 9.8 #gravitational acceleration [m/s2]

#######READ IN HYDRAULICS DATA FROM BRINKERHOFF ETAL 2019 TO GET SWOT VERSUS INEFFICIENT CHANNELS------------------------------
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
  filter(measured_rating_diff %in% c('Excellent', 'EXCL', 'GOOD', 'Good'))#, 'Fair', 'FAIR'))

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

#########QUICK DETOUR TO CALCULATE HG MODELS FOR FCO2 CALCULATIONS------------------------------
num <- nrow(data)
wid <- lm(log(width)~log(Qm3s), data=data)
dep <- lm(log(depth)~log(Qm3s), data=data)
models <- data.frame('model' = c('width', 'depth'), 'r2'=c(summary(wid)$r.squared, summary(dep)$r.squared), 'coef'=c(exp(wid$coefficient[1]), exp(dep$coefficient[1])), 'slope'=c(wid$coefficient[2], dep$coefficient[2]))
write.csv(models, 'cache\\AHG.csv')

#######CALCULATE CHANNEL GEOMETRY---------------------------------------------------------
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius, assuming rectangular channel [m]
data$ustar <- sqrt(g*data$Rh*data$slope) #friction/shear velocity [m/s]

#######SOME FLAGS FOR SWOT-OBSERVABLE RIVERS AND WHEN RH=H-------------------------------
data$flag_swot <- ifelse(data$width >= 100, 'SWOT', 'Small')

#########HYDRAULIC GEOMETRY OF SWOT-OBSERVABLE FLOWS
HG_swot <- round(mean((data[data$flag_swot == 'SWOT',]$Rh/data[data$flag_swot == 'SWOT',]$depth)),2)
n <- nrow(data[data$flag_swot == 'SWOT',])
df <- data.frame('mean'=HG_swot, 'n'=n)
write.csv(df, 'cache/k600_theory/HG_swot.csv')

#######LOG TRANSFORM SOME VARIABLES-----------------------------------------------------------
data$log_slope <- log(data$slope)

######## K600 MODELS--------------------
#hydraulically-wide rivers
data <- filter(data, study %in% c('Ulseth_etal_2019', 'Churchill_etal_1962', 'Owens_etal_1964'))
data$flag_hydraulicWide <- ifelse(data$Rh/data$depth >= 0.99, 'Rh=H', 'Rh=/=H')

hydraulicallyWide <- filter(data, flag_hydraulicWide == 'Rh=H')

#fit small-eddy model with surface manifest of bed dissipation
hydraulicallyWide$hydraulicallyWideModel <- g^(3/8)*hydraulicallyWide$slope^(3/8)*hydraulicallyWide$depth^(1/8)
lm_hydraulicallyWide_smallEddy_eS <- lm(k600~hydraulicallyWideModel+0, data=hydraulicallyWide)
hydraulicallyWide$k600_pred_wideHydraulics <- predict(lm_hydraulicallyWide_smallEddy_eS, hydraulicallyWide)

#plot model
plot_smallEddy_eS <- ggplot(hydraulicallyWide, aes(x=k600_pred_wideHydraulics, y=k600)) +
  geom_point(size=5, color='#bebada') +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_hydraulicallyWide_smallEddy_eS)$r.squared,2)), x = 1, y = 100, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(alpha*(gS)^{3/8}*H^{1/8}, ' [', m, '/', dy, ']'))),
       y = expression(bold(paste(k[600], ' [', m, '/', dy, ']'))))+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Small-eddy, ', epsilon==epsilon[S]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit reynolds model with surface manifest of bed dissipation
hydraulicallyWide$hydraulicallyWideModel <- (g*hydraulicallyWide$slope)^(9/16)*hydraulicallyWide$depth^(11/16)
lm_hydraulicallyWide_reynolds_eS <- lm(k600~hydraulicallyWideModel+0, data=hydraulicallyWide)
hydraulicallyWide$k600_pred_wideHydraulics <- (predict(lm_hydraulicallyWide_reynolds_eS, hydraulicallyWide))

#plot model
plot_reynolds_eS <- ggplot(hydraulicallyWide, aes(x=k600_pred_wideHydraulics, y=k600)) +
  geom_point(size=5, color='#bebada') +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_hydraulicallyWide_reynolds_eS)$r.squared,2)), x = 1, y = 100, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(beta*(gS)^{9/16}*H^{11/16}, ' [', m, '/', dy, ']'))),
       y = '')+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Reynolds extension, ', epsilon==epsilon[S]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit small-eddy model with form-drag dissipation
hydraulicallyWide$hydraulicallyWideModel <- (g*hydraulicallyWide$slope*hydraulicallyWide$Vms)^(1/4)
lm_hydraulicallyWide_smallEddy_eD <- lm(k600~hydraulicallyWideModel+0, data=hydraulicallyWide)
hydraulicallyWide$k600_pred_wideHydraulics <- predict(lm_hydraulicallyWide_smallEddy_eD, hydraulicallyWide)

#plot model
plot_smallEddy_eD <- ggplot(hydraulicallyWide, aes(x=k600_pred_wideHydraulics, y=k600)) +
  geom_point(size=5, color='#bebada') +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_hydraulicallyWide_smallEddy_eD)$r.squared,2)), x = 1, y = 100, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(alpha[1]*(gSU)^{1/4}, ' [', m, '/', dy, ']'))),
       y = expression(bold(paste(k[600], ' [', m, '/', dy, ']'))))+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Small-eddy, ', epsilon==epsilon[D]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#fit reynolds model with form-drag dissipation
hydraulicallyWide$hydraulicallyWideModel <- (g*hydraulicallyWide$slope)^(7/16)*hydraulicallyWide$Vms^(1/4)*hydraulicallyWide$depth^(9/16)
lm_hydraulicallyWide_reynolds_eD <- lm(k600~hydraulicallyWideModel+0, data=hydraulicallyWide)
hydraulicallyWide$k600_pred_wideHydraulics <- (predict(lm_hydraulicallyWide_reynolds_eD, hydraulicallyWide))

#plot model
plot_reynolds_eD <- ggplot(hydraulicallyWide, aes(x=k600_pred_wideHydraulics, y=k600)) +
  geom_point(size=5, color='#bebada') +
  geom_abline(linetype='dashed', color='darkgrey', size=1.5)+ #1:1 line
  annotate("text", label = paste0('r2: ', round(summary(lm_hydraulicallyWide_reynolds_eD)$r.squared,2)), x = 1, y = 100, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(beta[1]*(gS)^{7/16}*U^{1/4}*H^{9/16}, ' [', m, '/', dy, ']'))),
       y = '')+
  scale_y_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  scale_x_log10(limits=c(10^-1,10^2),
                breaks=c(0.1, 1, 10, 100),
                labels=c('0.1', '1', '10', '100'))+
  annotation_logticks()+
  ggtitle(expression(bold(paste('Reynolds extension, ', epsilon==epsilon[D]))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

k600_modelPlot <- plot_grid(plot_smallEddy_eS, plot_reynolds_eS, plot_smallEddy_eD, plot_reynolds_eD, ncol=2, label_size = 18, labels=c('a', 'b', 'c', 'd'))

ggsave('cache\\k600_theory\\k600Plot.jpg', k600_modelPlot, height=9, width=10)

#######WRITE reynolds MODEL TO FILE-------------------------------
models <- data.frame('name'=c('reynolds-eS', 'Small-eddy-eS', 'reynolds-eD', 'Small-eddy-eD'),
                     'r2'=c(summary(lm_hydraulicallyWide_reynolds_eS)$r.squared,
                            summary(lm_hydraulicallyWide_smallEddy_eS)$r.squared,
                            summary(lm_hydraulicallyWide_reynolds_eD)$r.squared,
                            summary(lm_hydraulicallyWide_smallEddy_eD)$r.squared),
                     'coef'=c(summary(lm_hydraulicallyWide_reynolds_eS)$coefficient[1],
                              summary(lm_hydraulicallyWide_smallEddy_eS)$coefficient[1],
                              summary(lm_hydraulicallyWide_reynolds_eD)$coefficient[1],
                              summary(lm_hydraulicallyWide_smallEddy_eD)$coefficient[1])
                     )
write.csv(models, 'cache\\k600_theory\\hydraulicWide_models.csv')


#######MONTE CARLO PROPOGATION OF UNCERTANTIES-------------------------
beta_sigma <- summary(lm_hydraulicallyWide_reynolds_eD)$sigma
u_sigma <- 0.3

set.seed(13)
n <- 10000 #MC sample size

output <- 1:nrow(hydraulicallyWide)
for (i in output) {
  log_knowns <- (7/16)*log(g) + (7/16)*log(hydraulicallyWide[i,]$slope) + (9/16)*log(hydraulicallyWide[i,]$depth)
  log_k600_pred <- rnorm(n, log(summary(lm_hydraulicallyWide_reynolds_eD)$coefficient[1]), log(beta_sigma)) + log_knowns + (1/4)*rnorm(n, log(hydraulicallyWide[i,]$Vms), u_sigma)
  output[i] <- sd(log_k600_pred)
}

mean(output)
hist(output)
