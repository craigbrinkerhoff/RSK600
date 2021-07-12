#################
##Description: Turbulence controls on k~turbulence dissipation relationship
##Creator: Craig Brinkerhoff
##Date: Sumeer 2021
#################

#######SETUP-----------------------------
library(tidyverse)
library(cowplot)
library(segmented)
theme_set(theme_cowplot())

setwd('C:\\Users\\craig\\Documents\\GitHub\\RSK600\\')

#some constants
g <- 9.8
cf <- 0.100 #Henderson found empirically that Cf=0.113, if we use the empirical Strickler relation for n (1923- different rivers), we get Cf=0.09. So, we used Cf=0.100. See the notes!
Sg <- 2.65 #specific gravity

#########READ IN HYDRAULICS DATA FROM BRINKERHOFF ETAL 2019------------------------------
data <- read.csv('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Published\\2019_AMHG_GRL\\working\\final\\field_measurements.csv')
nhd_med <- read.csv('C:\\Users\\craig\\Documents\\OneDrive - University of Massachusetts\\Published\\2019_AMHG_GRL\\working\\final\\NHD_join_table.csv')
data <- left_join(data, select(nhd_med, SLOPE, SOURCE_FEA), by = c("site_no" = "SOURCE_FEA"))

data <- filter(data, is.finite(chan_width) ==1) %>%
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
data$area <- data$chan_area * 0.092903 #ft2 to m2
data$width <- data$chan_width*0.305 #m
data$Vms <- data$chan_velocity*0.305 #m/s
data$depth <- data$chan_area / data$width
data$slope <- data$SLOPE


######READ IN DATA WITH K600-------------------
#Ulseth etal 2019
ulseth_data <- read.csv('data/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
ulseth_data <- ulseth_data[,-8]
ulseth_data$dataset <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$dataset == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$dataset))
ulseth_data <- filter(ulseth_data, is.na(width)==0) #filter out no width measurements
ulseth_data_s <- select(ulseth_data, 'width', 'Vms', 'depth', 'slope', 'k600')
ulseth_data_s$study <- 'Ulseth_etal_2019'

###############JOIN DATASETS---------------------------------
data <- select(data, 'width', 'Vms', 'depth', 'slope')
data$k600 <- NA
data$study <- 'Brinkerhoff_etal_2019'
data <- rbind(data, ulseth_data_s)

#########CALCULATE HYDRAULICS---------------------------------------------------------
data$eD <- g * data$slope * data$Vms #dissipation rate of surface turbulence originating from depth-scale form drag [W/kg]
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius [m]
data$ustar <- sqrt(g*data$Rh*data$slope)
data$eS <- (data$ustar^3)/(0.41*data$depth) #dissipation rate of surface turbulence originating from bed friction [W/kg]
data$Gtotal <- data$Vms*data$ustar^2

###########SOME FLAGS-------------------------------
data$flag_swot <- ifelse(data$width > 100, 'SWOT', 'Small')
data$flag_depth <- ifelse(round(data$Rh/data$depth, 3) >= 0.995, 'Rh=H', 'Rh=/=H')

#LOG TRANSFORM SOME VARIABLES-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_eS <- log(data$eS)
data$log_slope <- log(data$slope)

#####PARTIONING OF ED AND ES---------------------------------------------------------
#Similar to Brinkerhoff etal 2019 function. Calibrating a Shield's grain size using Manning's equation, Shield's parameter, and Henderson relation f~(D/depth)^(1/3)
#See math notes for the two velocity derivations and the math solving for the Henderson coefficient multiple ways. We ultimately use Cf=0.100

# #function to calibrate grain size
# grainSizer <- function(v,s, d,cf, Sg){
#   D <- 10^(seq(-5, 2, 0.005)) #virtual grain size [m] Tests 601 values
#   
#   #first method
#   d_prime <- v^(3/2) * (sqrt((8*g*s)/(cf)))^(-3/2) * (1/D^(1/6))^(-3/2) #[m] #eq 15
#   d_double_prime <- d - d_prime #[m]
#   inverse_shields_prime <- ((Sg-1)*D)/(s*d_prime) #eq 16
#   
#   C_1 <- ifelse(inverse_shields_prime < 1.6, 41, 29)#eqs 20-22
#   #ifelse(inverse_shields_prime < 17.8, 29, 23))
#   C_2 <- ifelse(inverse_shields_prime < 1.6, 1.26,0.5)
#   #ifelse(inverse_shields_prime < 17.8, 0.5, 0.415))
#   v_m_double_prime_1 <- C_1/inverse_shields_prime^C_2 #eq 19, approximation of Einstein & Barbarossa sediment regimes
#   
#   #2nd method
#   v_m_double_prime_2 <- v / sqrt(g*d_double_prime*s) #eq 14, closed form solution of v/m''
#   
#   #get calibrated value
#   diff <- abs(v_m_double_prime_1-v_m_double_prime_2) #minimal difference b/w skin friction depth estimates
#   D_fin <- D[which.min(diff)] #[m]
#   
#   #percent difference between estimates
#   perc <- diff/v_m_double_prime_2
#   perc <- min(perc, na.rm=T)*100
#   
#   return(c(D_fin, perc))
# }
# 
# #run grain sizer function
#grainResults <- mapply(grainSizer,data$Vms,data$slope, data$depth, cf, Sg)

# data$D35 <- grainResults[1,]
# data$grain_perc_diff <- grainResults[2,]
# summary(data$D35)
# summary(data$grain_perc_diff)
# data$log_D35 <- log(data$D35)
# 
# #calculate resistance partition coefficient using calibrated D35 (and its associated d')
# data$inverse_shields_prime <- data$Vms^(3/2) * (cf/(8*g*data$D35))^(3/4) * data$slope^(1/4) * (Sg-1)^(-1)
# data$d_prime <- data$inverse_shields_prime * data$slope^(-1) * (Sg-1) * data$D35 #calculated via shield's parameter and so D35. Remeber error is less than 5% so we good
# 
# data$f_prime <- cf *(data$D35/data$d_prime)^(1/3) #Henderson 1966
# 
# data$partition <- 1-((data$f_prime*data$Vms^2)/(8*g*data$Rh*data$slope)) #Moog & Jirka 19999
# data$partition <- ifelse(data$partition < 0, 0, data$partition)
# 
# #eM model vs k600------------------------------------------------
# data$eM <- (1-data$partition)*data$eS + data$partition*data$eD #combined Moog and Jirka 1999 macroroughness model
# data$log_eM <- log(data$eM)

##########GET ADDITIONAL HYDRAULICS FROM EM------------------------------
data$Td <- data$eD - data$Gtotal/data$depth #its impossible (at the free surface) for dissipation to be less than generation, so filter these out

##########PLOTS---------------
#TKE production vs TKE form-drag dissipation
#My analytical solution suggests these should be ~equal when depth equals Rh
#binning the Rh/H ratios for visualzing
data$figureFlag <- ifelse(round(data$Rh/data$depth, 3) == 1, 'Rh=H',
                          ifelse(round(data$Rh/data$depth, 3) >= 0.90, '90-99%',
                          ifelse(round(data$Rh/data$depth, 3) >= 0.80, '80-89%',
                                 ifelse(round(data$Rh/data$depth, 3) >= 0.70, '70-79%',
                                        ifelse(round(data$Rh/data$depth, 3) >= 0.60, '60-69%',
                                               ifelse(round(data$Rh/data$depth, 3) >= 0.50, '50-59%',
                                                      ifelse(round(data$Rh/data$depth, 3) >= 0.40, '40-49%',
                                                             ifelse(round(data$Rh/data$depth, 3) >= 0.30, '30-39%',
                                                                    ifelse(round(data$Rh/data$depth, 3) >=0.20, '20-29%',
                                                                           ifelse(round(data$Rh/data$depth, 3) >= 0.10, '10-19%', '0-9%'))))))))))
#create boxplots
boxes <- ggplot(data, aes(x=figureFlag, y=Td, fill=study)) +
  #geom_hline(yintercept = 1e-5, linetype='dashed', color='darkgrey', size=2)+
  geom_boxplot(size=1.2) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_brewer(palette='Accent', name='')+
  annotate('text', label='n=530,945', x='10-19%', y=10^-8, color='#7fc97f', size=8)+
  annotate('text', label='n=701', x='10-19%', y=10^-9, color='#beaed4', size=8)+
  ylab('eD-G [J/kg*s]') +
  xlab('Rh/H [%]') +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')
boxes

# scatter <- ggplot(data, aes(x=Rh/depth, y=Td, color=log10(Td))) +
#   geom_hline(yintercept = 1e-5, linetype='dashed', color='darkgrey', size=2)+
#   geom_point(size=5) +
#   geom_smooth(color='black', method='gam')+
#   annotate("text", label = '~0 J/kg*s', x = 0.6, y = 10^-4.75, size = 8, colour = "darkgrey")+
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   scale_color_distiller(palette='RdPu', name='log10 eD-G')+
#   ylab('e - G [J/kg*s]') +
#   xlab('Rh/H') +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=24,face="bold"),
#         legend.text = element_text(size=17),
#         legend.title = element_text(size=17, face='bold'))
# 
# legend <- get_legend(scatter)
# 
# finalPlot <- plot_grid(scatter + theme(legend.position='none'), boxes, ncol=2)
# finalPlot <- finalPlot + draw_grob(legend, x=0.125, y=-0.275, width=1, height=1)

ggplot(data, aes(x=width, linetype=flag_depth, color=study)) +
  stat_ecdf(size=2) +
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept = 100, linetype='dashed', size=2)+
  annotate('text', label='Observable\nvia SWOT', x=3000, y=0.3, size=8, color='black')+
  annotate('text', label='Not Observable\nvia SWOT', x=1, y=0.8, size=8, color='black')+
  scale_color_brewer(name='River Regime', palette = 'Accent')+
  scale_x_log10()


#k models plot----------------------------------------
data <- filter(data, study == 'Ulseth_etal_2019')
data$log_k600 <- log(data$k600)

#models
wideRegime <- data[data$flag_depth == 'Rh=H',]
lm_ustar <- lm(k600~ustar+0, data=wideRegime)
lm_eD <- lm(log10(k600)~log10(eD), data=wideRegime)
wideRegime$k600_pred_ustar <- predict(lm_ustar, wideRegime)
wideRegime$k600_pred_eD <- predict(lm_eD, wideRegime)

ustar_k600_wide <- ggplot(wideRegime) +
  geom_point(aes(x=ustar, y=k600, color=flag_swot), size=5) +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  scale_color_brewer(palette='Dark2')+
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.05, y = 17, size = 8, colour = "purple")+
  xlab('') +
  ylab('k600 [m/dy]') +
  ggtitle('U*-Wide')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

eD_k600_wide <- ggplot(wideRegime) +
  geom_point(aes(x=eD, y=k600, color=flag_swot), size=5) +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  scale_color_brewer(palette='Dark2')+
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-4, y = 17, size = 8, colour = "purple")+
  xlab('') +
  ylab('') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  ggtitle('eD-Wide')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

narrowRegime <- data[data$flag_depth == 'Rh=/=H',]
lm_ustar <- lm(k600~ustar+0, data=narrowRegime)
lm_eD <- lm(log10(k600)~log10(eD), data=narrowRegime)
narrowRegime$k600_pred_ustar <- predict(lm_ustar, narrowRegime)
narrowRegime$k600_pred_eD <- predict(lm_eD, narrowRegime)

ustar_k6002_narrow <- ggplot(narrowRegime) +
  geom_point(aes(x=ustar, y=k600, color=flag_swot), size=5) +
  geom_line(aes(x=ustar, y=k600_pred_ustar), size=2, color='black') +
  scale_color_brewer(palette='Dark2')+
  annotate("text", label = paste0('r2: ', round(summary(lm_ustar)$r.squared,2)), x = 0.1, y = 3000, size = 8, colour = "purple")+
  xlab('Ustar [m/s]') +
  ylab('k600 [m/dy]') +
  ggtitle('U*-Narrow')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

eD_k6002_narrow <- ggplot(narrowRegime) +
  geom_point(aes(x=eD, y=k600, color=flag_swot), size=5) +
  geom_line(aes(x=eD, y=10^k600_pred_eD), size=2, color='black') +
  scale_color_brewer(palette='Dark2', name='')+
  annotate("text", label = paste0('r2: ', round(summary(lm_eD)$r.squared,2)), x = 10^-4, y = 2000, size = 8, colour = "purple")+
  xlab('eD [J/kg*s]') +
  ylab('') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  ggtitle('eD-Narrow')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

plot_grid(ustar_k600_wide, eD_k600_wide, ustar_k6002_narrow, eD_k6002_narrow, ncol=2)

#save the actual model parameters
models_ustar <- group_by(data, flag_depth) %>%
  do(model=lm(k600~ustar+0, data=.)) %>%
  summarise(slope = model$coefficients[1],
            r2 = summary(model)$r.squared,
            log_SE = log(sqrt(deviance(model)/df.residual(model))),
            name = first(flag_depth))

models_eD <- group_by(data, flag_depth) %>%
  do(model=lm(log10(k600)~log10(eD), data=.)) %>%
  summarise(int = 10^model$coefficients[1],
            slope = model$coefficients[2],
            r2 = summary(model)$r.squared,
            log_SE = sqrt(deviance(model)/df.residual(model)),
            name = first(flag_depth))


#khat prior model
lm <- lm(log(ustar)~log_slope, data=wideRegime)
summary(lm)

#log_khat prior uncertainity
#Using error propogation for the ustar~slope model and the k600~ustar model
sqrt(sqrt(deviance(lm)/df.residual(lm))^2 + (models_ustar[2,]$SE)^2)

#plot
ggplot(wideRegime, aes(x=56.0294*exp(-0.69945)*slope^0.29587, y=56.0294*ustar)) +
  geom_point(size=3) +
  geom_abline()

