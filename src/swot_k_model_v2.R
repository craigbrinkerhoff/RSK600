#################
##Description: Developing gas exchange model for SWOT-observable rivers
##Creator: Craig Brinkerhoff
##Date: Fall 2021
#################

#######SETUP-----------------------------
library(tidyverse)
library(colorspace)
library(cowplot)
library(segmented)
theme_set(theme_cowplot())

setwd('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Ongoing Projects\\RSK600')

#some constants
g <- 9.8 #gravitational acceleration [m/s2]

#######READ IN HYDRAULICS DATA FROM BRINKERHOFF ETAL 2019 TO GET SWOT VERSUS INEFFICIENT CHANNELS------------------------------
data <- read.csv('data\\Brinkerhoff_etal_2019\\field_measurements.csv')
#nhd_med <- read.csv('data\\Brinkerhoff_etal_2019\\NHD_join_table.csv')
#data <- left_join(data, select(nhd_med, SLOPE, SOURCE_FEA), by = c("site_no" = "SOURCE_FEA"))

#some necessary filtering
data <- filter(data, is.finite(chan_width) ==1) %>%
  filter(is.finite(chan_velocity)==1) %>%
  filter(is.finite(chan_discharge)==1) %>%
  filter(is.finite(chan_area)==1) %>%
  filter(chan_width > 0) %>%
  filter(chan_velocity > 0) %>%
  filter(chan_discharge > 0) %>%
  filter(chan_area > 0) %>%
  filter(measured_rating_diff %in% c('Excellent', 'EXCL', 'GOOD', 'Good', 'Fair', 'FAIR'))# %>%
 # filter(is.finite(SLOPE)==1)

#imperial to metric
data$area <- data$chan_area * 0.092903 #ft2 to m2
data$width <- data$chan_width*0.305 #m
data$Vms <- data$chan_velocity*0.305 #m/s
data$depth <- data$area / data$width
data$slope <- NA #data$SLOPE
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

data <- filter(data, width <= bankful_width)
data <- filter(data, depth <= bankful_depth)

#get Dingman r channel shape (f / b in AHG parlance) and only keep measurements made at sites with robust r estimates > 1 (Dingman 2007)
# ahg <- group_by(data, site_no) %>%
#   do(model_w = lm(log(width)~log(Qm3s), data=.),
#      model_d = lm(log(depth)~log(Qm3s), data=.)) %>%
#   mutate(b = model_w$coefficient[2],
#          r2_w = summary(model_w)$r.squared,
#          f = model_d$coefficient[2],
#          r2_d = summary(model_d)$r.squared) %>%
#   select(site_no, b, f, r2_w, r2_d) %>%
#   mutate(r = f/b)
# data <- left_join(data, ahg, 'site_no')

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
data <- rbind(data, ulseth_data_s)

#######CALCULATE CHANNEL GEOMETRY---------------------------------------------------------
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius, assuming rectangular channel [m]
data$ustar <- sqrt(g*data$Rh*data$slope) #friction/shear velocity [m/s]

#######CALCULATE TURBULENCE TERMS-----------------------------------------
data$eD <- g * data$slope * data$Vms #dissipation rate of turbulence originating from depth-scale form drag [W/kg]
data$eS <- (data$ustar^3)/(data$depth) #dissipation rate of surface turbulence originating from bed shear [W/kg]

#######SOME FLAGS FOR SWOT-OBSERVABLE RIVERS AND WHEN RH=H-------------------------------
data$flag_swot <- ifelse(data$width >= 100, 'SWOT', 'Small')
data$flag_depth <- ifelse(round(data$Rh/data$depth, 3) >= 0.95, 'Rh=H', 'Rh=/=H') #set to 95% just to get idea of number of SWOT measurements that generally meet this condition
percs <- group_by(data, flag_swot, flag_depth) %>% summarise(n=n()) %>% group_by(flag_swot) %>% mutate(perc = 100 * n/sum(n))
write.csv(percs, 'cache/k600_theory/Rh_H_percents.csv')

#######LOG TRANSFORM SOME VARIABLES-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_eS <- log(data$eS)
data$log_slope <- log(data$slope)

####### K MODELS FITTING-------------------------------
data <- filter(data, study == 'Ulseth_etal_2019') #remove Brinkerhoff etal 2019 from dataset
calibrate_results <- data.frame('Rh_H'=NA, 'r2_ed_25'=NA, 'statistical_slope'=NA)

#loop through Rh/H thresholds and fit models to inefficient channels
for(i in seq(0.75, 0.999, 0.005)) {
  temp <- data #use temp storage of data
  temp$eD4 <- temp$eD^(1/4) #calculate Kolmogorov velocity scale
  temp$flag_depth <- ifelse(round(temp$Rh / temp$depth, 3) >= i, 'Rh=H', 'Rh=/=H') #calculate Rh_H threshold
  
  #fit models
  results <- group_by(temp, flag_depth) %>%
    do(model_ustar=lm(k600~ustar+0, data=.),
       model_eD4=lm(k600~eD4+0, data=.),
       model_eDgen=lm(log10(k600)~log10(eD), data=.)) %>%
    summarise(r2_ed_25 = summary(model_eD4)$r.squared,
              statistical_slope = model_eDgen$coefficients[2],
              name = first(flag_depth)) %>%
    filter(name == 'Rh=H') #only keep inefficient channel model
  output <- data.frame('Rh_H'=i, 'r2_ed_25'=results$r2_ed_25, 'statistical_slope'=results$statistical_slope)
  
  calibrate_results <- rbind(calibrate_results, output)
}
calibrate_results <- calibrate_results[-1,]

#plot results
forPlot <- gather(calibrate_results, key=key, value=value, c('r2_ed_25', 'statistical_slope'))
calibrationPlot <- ggplot(forPlot, aes(x=Rh_H, y=value, color=key)) +
  geom_point(size=3.5) +
  geom_line(size=1) +
  geom_hline(yintercept = 0.25, size=2, color='darkblue') +
  annotate("text", label = "r^2 * 'of ' * k[600] %prop% (epsilon[D])^(1/4)", x = 0.95, y = 0.07, size = 10, colour = "darkgreen", parse=TRUE)+
  annotate("text", label = "'b in ' * k[600] %prop% (epsilon[D])^b", x = 0.93, y = 0.70, size = 10, colour = "darkblue", parse=TRUE)+
  annotate("text", label = "Theoretical value (0.25)", x = 0.80, y = 0.30, size = 8, colour = "darkblue")+
  scale_color_manual(values=c('darkgreen', 'darkblue')) +
  ylim(0,1) +
  xlab('Inefficient channel threshold (minimum Rh/H)')+
  ylab('Value') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

####### K MODELS PLOTS-------------------------------
#create plots with the best fitting (and theoretically consistent) model
data$flag_depth <- ifelse(round(data$Rh/data$depth,3) >= calibrate_results[calibrate_results$r2_ed_25 == max(calibrate_results$r2_ed_25),]$Rh_H, 'Rh=H', 'Rh=/=H') #actually equal
data$log_k600 <- log(data$k600)
data$eD4 <- data$eD^(1/4) #calculate Kolmogorov velocity scale

#Inefficient river condition (Rh=H)-------------------------
inefficientRegime <- data[data$flag_depth == 'Rh=H',]
lm_ustar <- lm(k600~ustar+0, data=inefficientRegime) #depth-scale turbulence
lm_eD4 <- lm(k600~eD4+0, data=inefficientRegime) #Kolmogorov-scale turbulence
lm_eD <- lm(log10(k600)~log10(eD), data=inefficientRegime) #generalized Kolmorogorv-style model
inefficientRegime$k600_pred_ustar <- predict(lm_ustar, inefficientRegime)
inefficientRegime$k600_pred_eD4 <- predict(lm_eD4, inefficientRegime)
inefficientRegime$k600_pred_eD <- predict(lm_eD, inefficientRegime)

#depth-scale turbulence fitted model
ustar_k600_inefficient <- ggplot(inefficientRegime) +
  geom_point(aes(x=ustar, y=k600), size=5, color='#beaed4') +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.05, y = 17, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(U["*"], ' [m/s]'))),
       y = '')+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Kolmogorov-scale turbulence fitted model
eD4_k600_inefficient <- ggplot(inefficientRegime) +
  geom_point(aes(x=eD4, y=k600), size=5, color='darkgreen') +
  geom_line(aes(x=eD4, y=k600_pred_eD4), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD4)$r.squared,2)), x = 0.2, y = 17, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(epsilon[D]^(1/4), ' [', W^(1/4),'/',kg^(1/4), ']'))),
       y = '')+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#generalized Kolmorogorv-style fitted model: can we statistically recover 1/4?
generalized_k600_inefficient <- ggplot(inefficientRegime) +
  geom_point(aes(x=eD, y=k600), size=5, color='darkblue') +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-4, y = 17, size = 8, colour = "darkred")+
  labs(x = expression(bold(paste(epsilon[D], ' [W/kg]'))),
       y = expression(bold(paste(K[600], ' [m/dy]'))))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#efficient river condition: Rh =/= H NOT CURRENTLY BEING PLOTTED--------------------
efficientRegime <- data[data$flag_depth == 'Rh=/=H',]
lm_ustar <- lm(k600~ustar+0, data=efficientRegime) #depth-scale turbulence
lm_eD4 <- lm(k600~eD4+0, data=efficientRegime) #Kolmogorov-scale turbulence
lm_eD <- lm(log10(k600)~log10(eD), data=efficientRegime) #generalized Kolmorogorv-style model
efficientRegime$k600_pred_ustar <- predict(lm_ustar, efficientRegime)
efficientRegime$k600_pred_eD4 <- predict(lm_eD4, efficientRegime)
efficientRegime$k600_pred_eD <- predict(lm_eD, efficientRegime)

#depth-scale turbulence fitted model
ustar_k600_efficient <- ggplot(efficientRegime) +
  geom_point(aes(x=ustar, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.15, y = 2000, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(U["*"], ' [m/s]'))),
       y = expression(bold(paste(K[600], ' [m/dy]'))))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Kolmogorov-scale turbulence fitted model
eD4_k600_efficient <- ggplot(efficientRegime) +
  geom_point(aes(x=eD4, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=eD4, y=k600_pred_eD4), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD4)$r.squared,2)), x = 0.25, y = 2000, size = 8, colour = "purple")+
  labs(x = expression(bold(paste(epsilon[D]^(1/4), ' [', W^(1/4),'/',kg^(1/4), ']'))),
       y = '')+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#generalized Kolmorogorv-style fitted model: Can we statistically recover 1/4?
generalized_k600_efficient <- ggplot(efficientRegime) +
  geom_point(aes(x=eD, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-3.9, y = 100, size = 8, colour = "darkred")+
  labs(x = expression(bold(paste(epsilon[D], ' [W/kg]'))),
       y = '')+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#conceptual river efficiency models----------
efficientChannel_img <- ggdraw() +
  draw_image('cache/k600_theory/efficientChannel.jpg', scale=0.7)

inefficientChannel_img <- ggdraw() +
  draw_image('cache/k600_theory/inefficientChannel.jpg', scale=1)

arrow_img <- ggdraw() +
  draw_image('cache/k600_theory/arrow.jpg', scale=0.7)

#bring plots together & write to file
k600Plot <- plot_grid(efficientChannel_img, arrow_img, inefficientChannel_img, generalized_k600_inefficient, eD4_k600_inefficient, ustar_k600_inefficient, labels=c(NA, NA, NA, 'b', 'c', 'd'), ncol=3, label_size = 30)
k600Plot <- plot_grid(calibrationPlot, k600Plot, ncol=1, labels=c('a', NA), label_size = 30)
ggsave('cache\\k600_theory\\k600Plot.jpg', k600Plot, width=13, height=14)

#######WRITE THE ACTUAL MODEL PARAMETERS AND SCORES TO FILE-------------------------------
models_ustar <- group_by(data, flag_depth) %>%
  do(model=lm(k600~ustar+0, data=.)) %>%
  summarise(int = model$coefficients[1],
            slope = 1.0,
            r2 = summary(model)$r.squared,
            log_SE = log(sqrt(deviance(model)/df.residual(model))),
            name = first(flag_depth))

models_eD4 <- group_by(data, flag_depth) %>%
  do(model=lm(k600~eD4+0, data=.)) %>%
  summarise(int = model$coefficients[1],
            slope=0.25,
            r2 = summary(model)$r.squared,
            log_SE = sqrt(deviance(model)/df.residual(model)),
            name = first(flag_depth))

models_generalized <- group_by(data, flag_depth) %>%
  do(model=lm(log10(k600)~log10(eD), data=.)) %>%
  summarise(int = 10^model$coefficients[1],
            slope = model$coefficients[2],
            r2 = summary(model)$r.squared,
            log_SE = sqrt(deviance(model)/df.residual(model)),
            name = first(flag_depth))

write.csv(models_ustar, 'cache\\k600_theory\\ustar_models.csv')
write.csv(models_generalized, 'cache\\k600_theory\\ulseth_models.csv')
write.csv(models_eD4, 'cache\\k600_theory\\eD_models.csv')

#######PRINT MODELS FOR PRIOR SPECIFICATIONS. THESE ARE MANUALLY IMPLEMENTED WITHIN BIKER BUT ARE CALCULATED USING THIS DATASET--------------------------
#khat prior model for BIKER
lmPrior <- lm(log(ustar)~log_slope, data=inefficientRegime)
summary(lmPrior)

#log_khat prior uncertainity (Using error propogation and model standard errors for ustar~slope model and the k600~ustar model)
sqrt(sqrt(deviance(lmPrior)/df.residual(lmPrior))^2 + (models_ustar[2,]$log_SE)^2)

#Ulseth breakpoint model
lm <- lm(log_k600~log_eD, data=efficientRegime)
lm.seg <- segmented(lm, npsi=1)

summary(lm.seg)
