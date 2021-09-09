#################
##Description: Turbulence controls on k~turbulence dissipation relationship
##Creator: Craig Brinkerhoff
##Date: Summer 2021
#################

#######SETUP-----------------------------
library(tidyverse)
library(colorspace)
library(cowplot)
library(segmented)
theme_set(theme_cowplot())

setwd('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Ongoing Projects\\RSK600')

#some constants
g <- 9.8

#########READ IN HYDRAULICS DATA FROM BRINKERHOFF ETAL 2019------------------------------
data <- read.csv('data\\Brinkerhoff_etal_2019\\field_measurements.csv')
nhd_med <- read.csv('data\\Brinkerhoff_etal_2019\\NHD_join_table.csv')
data <- left_join(data, select(nhd_med, SLOPE, SOURCE_FEA), by = c("site_no" = "SOURCE_FEA"))

data <- filter(data, is.finite(chan_width) ==1) %>%
  filter(is.finite(chan_velocity)==1) %>%
  filter(is.finite(chan_discharge)==1) %>%
  filter(is.finite(chan_area)==1) %>%
  filter(chan_width > 0) %>%
  filter(chan_velocity > 0) %>%
  filter(chan_discharge > 0) %>%
  filter(chan_area > 0) %>%
  filter(measured_rating_diff %in% c('Excellent', 'EXCL', 'GOOD', 'Good', 'Fair', 'FAIR')) %>%
  filter(is.finite(SLOPE)==1)

#imperial to metric
data$area <- data$chan_area * 0.092903 #ft2 to m2
data$width <- data$chan_width*0.305 #m
data$Vms <- data$chan_velocity*0.305 #m/s
data$depth <- data$chan_area / data$width
data$slope <- data$SLOPE
data$Qm3s <- data$chan_discharge * 0.0283 #ft3/s to m3/s

######READ IN FIELD DATA WITH K600-------------------
#Ulseth etal 2019
ulseth_data <- read.csv('data/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM")
ulseth_data <- ulseth_data[,-8]
ulseth_data$dataset <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$dataset == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$dataset)) #contains data from 5 studies
ulseth_data <- filter(ulseth_data, is.na(width)==0) #filter out no width measurements
ulseth_data_s <- select(ulseth_data, 'width', 'Vms', 'depth', 'slope', 'k600', 'Qm3s')
ulseth_data_s$study <- 'Ulseth_etal_2019'

###############JOIN DATASETS---------------------------------
data <- select(data, 'width', 'Vms', 'depth', 'slope', 'Qm3s')
data$k600 <- NA
data$study <- 'Brinkerhoff_etal_2019'
data <- rbind(data, ulseth_data_s)

#########CALCULATE HYDRAULICS---------------------------------------------------------
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius [m]
data$ustar <- sqrt(g*data$Rh*data$slope) #friction/shear velocity [m/s]

#########CALCULATE TURBULENCE TERMS-----------------------------------------
data$eD <- g * data$slope * data$Vms #dissipation rate of surface turbulence originating from depth-scale form drag [W/kg]
data$eS <- (data$ustar^3)/(data$depth) #dissipation rate of surface turbulence originating from bed friction [W/kg]
data$Gtotal <- (data$Vms*data$ustar^2)/data$depth #depth-scale TKE produced, following Nakagawa & Nexu (1993) [W/kg]
data$Td <- data$eD - data$Gtotal #Turbulent diffusion of TKE from bed to the free surface, solved by closing the energy balance [W/kg]

###########SOME FLAGS FOR SWOT-OBSERVABLE RIVERS AND WHEN RH=H-------------------------------
data$flag_swot <- ifelse(data$width >= 100, 'SWOT', 'Small')
data$flag_depth <- ifelse(round(data$Rh/data$depth, 3) >= 0.995, 'Rh=H', 'Rh=/=H')

#LOG TRANSFORM SOME VARIABLES-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_eS <- log(data$eS)
data$log_slope <- log(data$slope)

##########TKE PLOTS----------------------------------
#binning the Rh/H ratios and Q for ease in visualzing
data$figureFlag <- ifelse(round(data$Rh/data$depth, 3) >= 0.995, 'Rh=H',
                          ifelse(round(data$Rh/data$depth, 3) >= 0.90, '90-99',
                            ifelse(round(data$Rh/data$depth, 3) >= 0.80, '80-89',
                                 ifelse(round(data$Rh/data$depth, 3) >= 0.70, '70-79',
                                        ifelse(round(data$Rh/data$depth, 3) >= 0.60, '60-69',
                                               ifelse(round(data$Rh/data$depth, 3) >= 0.50, '50-59',
                                                      ifelse(round(data$Rh/data$depth, 3) >= 0.40, '40-49',
                                                             ifelse(round(data$Rh/data$depth, 3) >= 0.30, '30-39',
                                                                    ifelse(round(data$Rh/data$depth, 3) >=0.20, '20-29',
                                                                           ifelse(round(data$Rh/data$depth, 3) >= 0.10, '10-19', '0-9'))))))))))
data$QFlag <- ifelse(data$Qm3s <= 0.0001, '0.0001',
                     ifelse(data$Qm3s <= 0.001, '0.001',
                            ifelse(data$Qm3s <= 0.01, '0.01',
                                   ifelse(data$Qm3s <= 0.1, '0.1',
                                          ifelse(data$Qm3s <= 1, '1',
                                                 ifelse(data$Qm3s <= 10, '10',
                                                        ifelse(data$Qm3s <= 100, '100',
                                                               ifelse(data$Qm3s <= 1000, '1000', '1000+'))))))))

#useful stats across both datasets tested
turbulence_stats <- group_by(data, figureFlag) %>% 
  summarise(meanTd = mean(Td),
            medianTd = median(Td),
            sdTd = sd(Td),
            n=n())
write.csv(turbulence_stats, 'cache\\k600_theory\\turbulence_stats.csv')

#plot TKE G/eD discrepancy versus channel shape
boxes_TKE <- ggplot(data, aes(x=flag_depth, y=Td, fill=flag_depth)) +
  geom_boxplot(size=1.2) +
  scale_fill_brewer(palette='Accent', name='')+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotate('text', label='n = 475,142', x='Rh=H', y=10^0, color='darkblue', size=10)+
  ylab('eD-G [J/kg*s]') +
  xlab('') +
  theme(axis.text.y=element_text(size=19),
        axis.text.x = element_text(size=24, face='bold'),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none') #c(0.55, 0.15)

#Q versus channel shape plot
w_d_stats <- group_by(data, QFlag) %>% 
  summarise(meanWD = mean(width/depth),
            sdWD = sd(width/depth),
            n=n())

wdPlot <- ggplot(data, aes(y=width/depth, x=QFlag)) +
  geom_boxplot(size=1.2, fill='lightblue') +
  geom_hline(yintercept = 1, linetype='dashed', size=1.2) +
  coord_cartesian(ylim=c(10^-2, 10^2))+
  annotate('text', label='(Axes cropped \nfor visualization)', x='0.01', y=10^-1.8, color='darkblue', size=8)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylab('Width/Depth') +
  xlab('Q [m3/s]') +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Td versus versus depth in hypotetical water-column (via Branch etal 2021 and Nexu 1993)
df <- data.frame('z_H'=seq(0,1,0.01))
df$tke_5 <- 4.78*sqrt(9.8*0.0001*5)^2*exp(-2*df$z_H) #5m depth
df$tke_3 <- 4.78*sqrt(9.8*0.0001*3)^2*exp(-2*df$z_H) #3m depth
df$tke_1 <- 4.78*sqrt(9.8*0.0001*1)^2*exp(-2*df$z_H) #1m depth
df$tke_05 <- 4.78*sqrt(9.8*0.0001*0.5)^2*exp(-2*df$z_H) #1/2m depth
df <- gather(df, key=key, value=value, c('tke_05', 'tke_1', 'tke_3', 'tke_5'))

#get tke gradient
df <- group_by(df, key) %>%
  mutate(dTKE = max(value)-min(value))

tke_depthPlot <- ggplot(df, aes(x=value, y=z_H, color=key, size=dTKE)) +
  geom_line() +
  geom_hline(yintercept = 1, size=3)+
  annotate("text", label = 'Free Surface', x = 0.0175, y = 0.97, size = 7)+
  geom_hline(yintercept = 0, size=3)+
  annotate("text", label = 'Riverbed', x = 0.0175, y = 0.03, size = 7)+
  scale_color_discrete_sequential(palette = "Purples 3", nmax = 9, order = 6:9, name='Depth of River Flow [m]', labels=c('0.5m', '1m', '3m', '5m'))+
  scale_size(name='TKE Gradient [J/kg]')+
  xlab('TKE [J/kg]') +
  ylab('Distance above river bed') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=17, face='bold'),
        legend.position = c(0.55, 0.65)) +
  guides(color = guide_legend(override.aes = list(size=5)))

#bring it alllll together and plot
plots_2 <- plot_grid(wdPlot, tke_depthPlot, ncol=2, labels=c('b','c'), label_size = 20)
plot_1 <- plot_grid(boxes_TKE, plots_2, ncol=1, labels=c('a',NA), label_size = 20)
ggsave('cache\\k600_theory\\turbulence.jpg', plot_1, width=14, height=15)

#confirm that depth-scale and Kolmogorov-scale turbulence is ~similar when Rh=H (used later in the k600 plot but need in memory now)--------------
ratio_stats <- group_by(data, figureFlag) %>% 
  summarise(meanratio = mean(((1e-6*eD)^0.25)/ustar),
            medianRatio = median(((1e-6*eD)^0.25)/ustar),
            n=n())
write.csv(ratio_stats, 'cache\\k600_theory\\V_ratio_stats.csv')

velocityScalePlot <- ggplot(data, aes(x=flag_depth, y=(((1e-6*eS)^0.25)/ustar), fill=flag_depth)) + #kinematic viscosity is assumed 1e-6 m2/s because no water temperature data is available
  geom_boxplot(size=1.2) +
  scale_fill_brewer(palette='Accent', name='')+
  geom_hline(yintercept = 1, linetype='dashed', size=2) +
  annotate('text', label='n = 475,142', x='Rh=H', y=10^-2, color='darkblue', size=10)+
  annotate('text', parse=TRUE, label = "frac((nu*epsilon[D])^{1/4}, Ustar)", x='Rh=/=H', y=10^0.5, color='darkblue', size=10)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylab('Characteristic turbulent \nvelocity ratio [%]') +
  xlab('') +
  theme(axis.text.y=element_text(size=19),
        axis.text.x = element_text(size=24, face='bold'),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Plot distributions of Rh/H ratios for SWOT-observable and non-SWOT observable rivers
swotPlot <- ggplot(data, aes(x=flag_depth, y=Qm3s, fill=flag_depth)) + #kinematic viscosity is assumed 1e-6 m2/s because no water temperature data is available
  geom_boxplot(size=1.2) +
  scale_fill_brewer(palette='Accent', name='')+
 # geom_hline(yintercept = 50, linetype='dashed', size=2) +
#  geom_hline(yintercept = 100, linetype='dashed', size=2) +
  annotate('text', label='n = 475,142', x='Rh=H', y=10^-2, color='darkblue', size=10)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylab('River Width [m]') +
  xlab('') +
  theme(axis.text.y=element_text(size=19),
        axis.text.x = element_text(size=24, face='bold'),
        axis.title=element_text(size=24,face="bold"),
        legend.position = 'none')
ggsave('cache\\k600_theory\\swotPlot.jpg', swotPlot, width=8, height=8)

###################### K MODELS PLOTS-------------------------------
data <- filter(data, study == 'Ulseth_etal_2019') #remove Brinkerhoff etal 2019 from dataset
data$log_k600 <- log(data$k600)
data$eD4 <- data$eD^(1/4) #calculate Kolmogorov velocity scale

#Large river model (Rh=H)
wideRegime <- data[data$flag_depth == 'Rh=H',]
lm_ustar <- lm(k600~ustar+0, data=wideRegime) #depth-scale turbulence
lm_eD4 <- lm(k600~eD4+0, data=wideRegime) #Kolmogorov-scale turbulence
lm_eD <- lm(log10(k600)~log10(eD), data=wideRegime) #generalized Kolmorogorv-style model
wideRegime$k600_pred_ustar <- predict(lm_ustar, wideRegime)
wideRegime$k600_pred_eD4 <- predict(lm_eD4, wideRegime)
wideRegime$k600_pred_eD <- predict(lm_eD, wideRegime)

#plot model results
#depth-scale turbulence
ustar_k600_wide <- ggplot(wideRegime) +
  geom_point(aes(x=ustar, y=k600), size=5, color='#beaed4') +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.05, y = 17, size = 8, colour = "purple")+
  xlab('') +
  ylab('k600 [m/dy]') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Kolmogorov-scale turbulence
eD4_k600_wide <- ggplot(wideRegime) +
  geom_point(aes(x=eD4, y=k600), size=5, color='#beaed4') +
  geom_line(aes(x=eD4, y=k600_pred_eD4), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD4)$r.squared,2)), x = 0.2, y = 17, size = 8, colour = "purple")+
  xlab('') +
  ylab('') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#generalized Kolmorogorv-style model
ulseth_k600_wide <- ggplot(wideRegime) +
  geom_point(aes(x=eD, y=k600), size=5, color='#beaed4') +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-4, y = 17, size = 8, colour = "darkred")+
  xlab('') +
  ylab('') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
               labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
               labels = scales::label_comma(drop0trailing = TRUE))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#narrow river condition: Rh =/= H
narrowRegime <- data[data$flag_depth == 'Rh=/=H',]
lm_ustar <- lm(k600~ustar+0, data=narrowRegime) #depth-scale turbulence
lm_eD4 <- lm(k600~eD4+0, data=narrowRegime) #Kolmogorov-scale turbulence
lm_eD <- lm(log10(k600)~log10(eD), data=narrowRegime) #generalized Kolmorogorv-style model
narrowRegime$k600_pred_ustar <- predict(lm_ustar, narrowRegime)
narrowRegime$k600_pred_eD4 <- predict(lm_eD4, narrowRegime)
narrowRegime$k600_pred_eD <- predict(lm_eD, narrowRegime)

#depth-scale turbulence
ustar_k600_narrow <- ggplot(narrowRegime) +
  geom_point(aes(x=ustar, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.15, y = 2000, size = 8, colour = "purple")+
  xlab('Ustar [m/s]') +
  ylab('k600 [m/dy]') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#Kolmogorov-scale turbulence
eD4_k600_narrow <- ggplot(narrowRegime) +
  geom_point(aes(x=eD4, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=eD4, y=k600_pred_eD), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD4)$r.squared,2)), x = 0.25, y = 2000, size = 8, colour = "purple")+
  xlab('eD^(1/4) [J/kg*s]^(1/4)') +
  ylab('') +
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#generalized Kolmorogorv-style model
ulseth_k600_narrow <- ggplot(narrowRegime) +
  geom_point(aes(x=eD, y=k600), size=5, color='#7fc97f') +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-4, y = 100, size = 8, colour = "darkred")+
  xlab('eD [J/kg*s]') +
  ylab('') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::label_comma(drop0trailing = TRUE))+
  theme(axis.text=element_text(size=19),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

#bringin' it allllll together
#text1 <- ggdraw() + draw_label("Big River Condition", fontface = 'bold', x = 0.5, y=0.65, size = 20)
#text2 <- ggdraw() + draw_label("All Other Rivers", fontface = 'bold', x = 0.5, y=0.75, size=20)

#save to file
k600Plot <- plot_grid(ustar_k600_wide, eD4_k600_wide, ulseth_k600_wide, ustar_k600_narrow, eD4_k600_narrow, ulseth_k600_narrow, labels=c('a', 'b', 'c', 'd', 'e', 'f'), ncol=3, label_size = 20)
k600Plot <- plot_grid(k600Plot, velocityScalePlot, ncol=1, labels=c(NA, 'g'), label_size = 20)
ggsave('cache\\k600_theory\\k600Plot.jpg', k600Plot, width=13, height=13)

###########SAVE THE ACTUAL MODELS-------------------------------
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

models_ulseth <- group_by(data, flag_depth) %>%
  do(model=lm(log10(k600)~log10(eD), data=.)) %>%
  summarise(int = 10^model$coefficients[1],
            slope = model$coefficients[2],
            r2 = summary(model)$r.squared,
            log_SE = sqrt(deviance(model)/df.residual(model)),
            name = first(flag_depth))

write.csv(models_ustar, 'cache\\k600_theory\\ustar_models.csv')
write.csv(models_ulseth, 'cache\\k600_theory\\ulseth_models.csv')
write.csv(models_eD4, 'cache\\k600_theory\\eD_models.csv')

############PRINT MODELS FOR PRIOR SPECIFICATIONS. THESE ARE MANUALLY IMPLEMENTED WITHIN BIKER BUT ARE CALCULATED USING THIS DATASET--------------------------
#khat prior model for BIKER
lmPrior <- lm(log(ustar)~log_slope, data=wideRegime)
summary(lmPrior)

#log_khat prior uncertainity (Using error propogation and model standard errors for ustar~slope model and the k600~ustar model)
sqrt(sqrt(deviance(lmPrior)/df.residual(lmPrior))^2 + (models_ustar[2,]$log_SE)^2)

#Ulseth breakpoint model
lm <- lm(log_k600~log_eD, data=narrowRegime)
lm.seg <- segmented(lm, npsi=1)

summary(lm.seg)


########
##SCRAPS------------------
##########

#Calculate hydraulic geometry
# dfHG <- group_by(data, site_no) %>%
#   filter(n() >= 20) %>%
#   do(model_w = lm(log(width)~log(discharge), data=.),
#      model_d = lm(log(depth)~log(discharge), data=.)) %>%
#   summarise(site_no = mean(site_no),
#             b = model_w$coefficient[2],
#             f = model_d$coefficient[2],
#             r2_width = summary(model_w)$r.squared,
#             r2_depth = summary(model_d)$r.squared) %>%
#   filter(r2_width > 0.6 & r2_depth > 0.6)
# dfHG$r <- round(dfHG$f / dfHG$b, 1)
# 
# data <- left_join(data, dfHG, by='site_no')
# data <- filter(data, is.na(r) == 0 & r > 1)
# 
# data$Hm <- ((data$r+1)/data$r)*data$depth #cross-section max depth per Dingman 2007
# data$Rh <- (data$width*data$depth)/(data$width + (8/3)*(data$Hm^2/data$width)) #Obtained using Dingman model for cross-section max depth (accountsfor channel shape)


#####PARTIONING OF ED AND ES---------------------------------------------------------
# cf <- 0.100 #Henderson found empirically that Cf=0.113, if we use the empirical Strickler relation for n (1923- different rivers), we get Cf=0.09. So, we used Cf=0.100. See the notes!
# Sg <- 2.65 #specific gravity
# 
# #Similar to Brinkerhoff etal 2019 function. Calibrating a Shield's grain size using Manning's equation, Shield's parameter, and Henderson relation f~(D/depth)^(1/3)
# #See math notes for the two velocity derivations and the math solving for the Henderson coefficient multiple ways. We ultimately use Cf=0.100
# 
#  #function to calibrate grain size
#  grainSizer <- function(v,s, d,cf, Sg){
#    D <- 10^(seq(-5, 2, 0.005)) #virtual grain size [m] Tests 601 values
# 
#    #first method
#    d_prime <- v^(3/2) * (sqrt((8*g*s)/(cf)))^(-3/2) * (1/D^(1/6))^(-3/2) #[m] #eq 15
#    d_double_prime <- d - d_prime #[m]
#    inverse_shields_prime <- ((Sg-1)*D)/(s*d_prime) #eq 16
# 
#    C_1 <- ifelse(inverse_shields_prime < 1.6, 41, 29)#eqs 20-22
#    #ifelse(inverse_shields_prime < 17.8, 29, 23))
#    C_2 <- ifelse(inverse_shields_prime < 1.6, 1.26,0.5)
#    #ifelse(inverse_shields_prime < 17.8, 0.5, 0.415))
#    v_m_double_prime_1 <- C_1/inverse_shields_prime^C_2 #eq 19, approximation of Einstein & Barbarossa sediment regimes
# 
#    #2nd method
#    v_m_double_prime_2 <- v / sqrt(g*d_double_prime*s) #eq 14, closed form solution of v/m''
# 
#    #get calibrated value
#    diff <- abs(v_m_double_prime_1-v_m_double_prime_2) #minimal difference b/w skin friction depth estimates
#    D_fin <- D[which.min(diff)] #[m]
# 
#    #percent difference between estimates
#    perc <- diff/v_m_double_prime_2
#    perc <- min(perc, na.rm=T)*100
# 
#    return(c(D_fin, perc))
#  }
# 
#  #run grain sizer function
# grainResults <- mapply(grainSizer,data$Vms,data$slope, data$depth, cf, Sg)
# 
#  data$D35 <- grainResults[1,]
#  data$grain_perc_diff <- grainResults[2,]
#  summary(data$D35)
#  summary(data$grain_perc_diff)
#  data$log_D35 <- log(data$D35)
#  

#calculate resistance partition coefficient using calibrated D35 (and its associated d')
 # data$inverse_shields_prime <- data$Vms^(3/2) * (cf/(8*g*data$D35))^(3/4) * data$slope^(1/4) * (Sg-1)^(-1)
 # data$d_prime <- data$inverse_shields_prime * data$slope^(-1) * (Sg-1) * data$D35 #calculated via shield's parameter and so D35. Remeber error is less than 5% so we good
 # 
 # data$f_prime <- cf *(data$D35/data$d_prime)^(1/3) #Henderson 1966
 # 
 # data$partition <- 1-((data$f_prime*data$Vms^2)/(8*g*data$Rh*data$slope)) #Moog & Jirka 19999
 # data$partition <- ifelse(data$partition < 0, 0, data$partition)
 # 
