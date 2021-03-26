###########################
#Description: Modeling k600 from hydraulics/roughness. Using open-channel flow to explain k600 scaling dyanmics with surface turbulence
#Creator: Craig Brinkerhoff
#Spring 2021
###########################

#User inputs------------------
cf <- 0.100 #Henderson found empirically that Cf=0.113, if we use the empirical Strickler relation for n (1923- different rivers), we get Cf=0.09. So, we used Cf=0.100. See the notes!
Sg <- 2.65 #specific gravity

#load model settings and packages-------------------------------
source(here::here('scripts' , 'inputs.R'))

#Ulseth etal 2019------------------------------------------------------------
ulseth_data <- read.csv('inputs//k600//Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
ulseth_data <- ulseth_data[,-8]
ulseth_data$data <- ifelse(ulseth_data$data == 'This Study', 'Ulseth et al. 2019', as.character(ulseth_data$data)) #fix some labeling
ulseth_data$data <- ifelse(ulseth_data$data == 'Raymond et al. 2012', 'Raymond et al. 2012 [5 studies]', as.character(ulseth_data$data))

data <- ulseth_data
data <- filter(data, is.na(width)==0) #filter out no width measurements

#calculate eD and other hydraulics---------------------------------------------------------
data$eD <- g * data$slope * data$Vms #dissipation rate of surface turbulence originating from depth-scale form drag [W/kg]
data$eS <- (sqrt(data$depth*g*data$slope))^3/data$depth #dissipation rate of surface turbulence originating from bed friction [W/kg]
data$Rh <- (data$depth*data$width)/(data$width + 2*data$depth) #hydraulic radius [m]
data$area <- data$width*data$depth #channel area [m2]

#log transform some variables-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_k600 <- log(data$k600)
data$log_slope <- log(data$slope)
data$log_area <- log(data$area)


#partitoning of eD and eS---------------------------------------------------------
#Similar to Brinkerhoffetal 2019 function. Calibrating a Shield's grain size using Manning's equation, Shield's parameter, and Henderson relation f~(D/depth)^(1/3)
  #See math notes for the two velocity derivations and the math solving for the Henderson coefficient multiple ways. We ultimately use Cf=0.100
grainSizer <- function(v,s, cf, Sg){
  D <- 10^(seq(-10, 2, 0.02)) #virtual grain size [m] Tests 601 values

  #test grain sizes through two methods to solve for d'
  #ASSUMPTIONS: Manning's parameters + Henderson relation
  d_prime_1 <- v^(3/2) * (sqrt((8*g*s)/(cf)))^(-3/2) * (1/D^(1/6))^(-3/2)

  #ASSUMPTIONS: Manning's parameters + Henderson relation + Shield's parameter for sediment entrainment (shields is inverted for algebraic convenience)
  inverse_shields_prime <- (1/v)*sqrt((8*g*D)/cf)*(((Sg-1)^(2/3))/(s^(1/6)))
  d_prime_2 <- inverse_shields_prime^(-1) * s^(-1) * ((Sg-1)*D)

  #get calibrated value
  diff <- abs(d_prime_1-d_prime_2) #minimal difference b/w skin friction depth estimates
  D_fin <- D[which.min(diff)]

  return(D_fin)
}

#run grain sizer
data$D <- mapply(grainSizer,data$Vms,data$slope, cf, Sg)
summary(data$D)
data$log_D <- log(data$D)

#calculate resistance partition coefficient using calibrated D and d'
data$inverse_shields_prime <- data$Vms^(3/2) * (cf/(8*g*data$D))^(3/4) * data$slope^(1/4) * (Sg-1)^(-1) #in the actual implementation we use all three assumptions
data$d_prime <- data$inverse_shields_prime * data$slope^(-1) * (Sg-1) * data$D #ASSUMPTIONS: Manning's parameters + Henderson relation + Shield's parameter for sediment entrainment

data$f_prime <- cf *(data$D/data$d_prime)^(1/3) #Henderson 1966

data$partition <- 1-((data$f_prime*data$Vms^2)/(8*g*data$Rh*data$slope)) #Moog & Jirka 19999
data$partition <- ifelse(data$partition < 0, 0, data$partition)

#eM model vs k600------------------------------------------------
data$eM <- (1-data$partition)*data$eS + data$partition*data$eD #combined Moog and Jirka macroroughness model
data$log_eM <- log(data$eM)

#predict k600 using some models--------------------------------------
#3 piece eM fit
linearModel_eM <- lm(log_k600~log_eM, data=data)
lm_craig_eM <- segmented(linearModel_eM, npsi = 2) # two breakpoints
r2_craig_eM <- round(summary(lm_craig_eM)$r.squared, 2)
summary(lm_craig_eM)
data$k600_eM <- exp(predict(lm_craig_eM, data))

#Craig's manual breakpoints
data$eM_regime_craig <- ifelse(data$eM < 10^-3.3, 'Small',
                          ifelse(data$eM < 10^-1, 'Medium', 'High'))

#3 piece D fit
linearModel_D <- lm(log_k600~log_D, data=data)
lm_craig_D <- segmented(linearModel_D, npsi = 1) # two breakpoints
r2_craig_D <- round(summary(lm_craig_D)$r.squared, 2)
summary(lm_craig_D)
data$k600_D <- exp(predict(lm_craig_D, data))

#3 piece area fits
linearModel_area <- lm(log_area~log_eM, data=data)
lm_craig_area <- segmented(linearModel_area, npsi = 2) # two breakpoints
r2_craig_area <- round(summary(lm_craig_area)$r.squared, 2)
data$area_pred <- exp(predict(lm_craig_area, data))

#2 piece slope fits
linearModel_slope <- lm(log_slope~log_eM, data=data)
lm_craig_slope <- segmented(linearModel_slope, npsi = 1) # two breakpoints
r2_craig_slope <- round(summary(lm_craig_slope)$r.squared, 2)
data$slope_pred <- exp(predict(lm_craig_slope, data))

#plots-------------------------------------------------------------------
eD_eS <- ggplot(data, aes(x=eD, y=eS, color=slope)) +
  geom_point(size=4, alpha=0.50) +
  geom_abline(linetype='dashed', size=2, color='black')+
  annotate("text", label = "eD ~should be \ngreater than eS", x = 10^-4, y = 10^0, size = 5, colour = "black") +
  xlab('eD [W/kg]') +
  ylab('eS [W/kg]') +
  scale_color_distiller(palette='Spectral')+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//roughness_theory//eD_eS.jpg', eD_eS, width=9, height=7)

eM_k600 <- ggplot(data, aes(x=eM, y=k600)) +
  geom_point(size=4, alpha=0.50) +
  geom_line(aes(x=eM, y=k600_eM), color='darkgreen', size=3) +
  geom_hline(yintercept=35, color='black', linetype='dashed', size=2)+
  annotate("text", label = "Ulseth: ~maximum \ndiffusive gas exchange (35 m/dy)", x = 10^-3, y = 10^2, size = 5, colour = "black") +
  xlab('eM [W/kg]') +
  ylab('k600 [m/dy]') +
  geom_richtext(aes(x=10^-4, y=10^3), label=paste0('r<sup>2</sup>: ', r2_craig_eM), color='darkgreen') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'bottom')
ggsave('outputs//k600//roughness_theory//eM_k600.jpg', eM_k600, width=9, height=7)

eM_partition <- ggplot(data, aes(x=eM, y=partition*100, color=eM_regime_craig)) +
  geom_point(size=4, alpha=0.50) +
#  annotate("text", label = "Regression Breakpoints", x = 10^-4, y = 90, size = 5, colour = "black") +
#  geom_vline(xintercept=0.00205, linetype='dashed', color='black', size=2) +
#  geom_vline(xintercept=0.1126, linetype='dashed', color='black', size=2) +
  xlab('eM [W/kg]') +
  ylab('Turbulence due to \nDepth-Scale Form Drag [%]') +
  scale_color_brewer(palette='Set2')+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//roughness_theory//eM_partition.jpg', eM_partition, width=9, height=7)

eM_D <- ggplot(data, aes(x=D*1000, y=k600, color=slope)) +
  geom_point(size=4, alpha=0.75) +
  geom_line(aes(x=D*1000, y=k600_D), color='darkgreen', size=3) +
  geom_hline(yintercept=35, color='black', linetype='dashed', size=2)+
  annotate("text", label = "Ulseth: ~maximum \ndiffusive gas exchange (35 m/dy)", x = 10^-3, y = 10^2, size = 5, colour = "black") +
  xlab('D [mm]') +
  ylab('k600 [m/dy]') +
  geom_richtext(aes(x=10^-4, y=10^3), label=paste0('r<sup>2</sup>: ', r2_craig_D), color='darkgreen') +
  scale_color_distiller(palette='Spectral')+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//roughness_theory//eM_D.jpg', eM_D, width=9, height=7)

partitionPlot <- ggplot(data, aes(x=partition*100, color=eM_regime_craig)) +
    stat_ecdf(size=2) +
    xlab('Turbulence due to \nDepth-Scale Form Drag [%]') +
    scale_color_brewer(palette='Set2', name='eM breakpoint regime')+
    theme(axis.text=element_text(size=20),
      axis.title=element_text(size=24,face="bold"),
      legend.text = element_text(size=17),
      legend.title = element_text(size=17, face='bold'),
     legend.position='none')
ggsave('outputs//k600//roughness_theory//partitionPlot.jpg', partitionPlot, width=9, height=7)

eM_area <- ggplot(data, aes(x=eM, y=area)) +
  geom_point(size=4, alpha=0.50) +
  geom_line(aes(x=eM, y=area_pred), color='darkgreen', size=3) +
  xlab('eM [W/kg]') +
  ylab('Channel Area [m2]') +
  geom_richtext(aes(x=10^-4, y=10^0), label=paste0('r<sup>2</sup>: ', r2_craig_area), color='darkgreen') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

eM_slope <- ggplot(data, aes(x=eM, y=slope)) +
  geom_point(size=4, alpha=0.50) +
  geom_line(aes(x=eM, y=slope_pred), color='darkgreen', size=3) +
  xlab('eM [W/kg]') +
  ylab('Channel Slope') +
  geom_richtext(aes(x=10^-4, y=10^0), label=paste0('r<sup>2</sup>: ', r2_craig_slope), color='darkgreen') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'),
        legend.position = 'none')

eM_channel <- plot_grid(eM_area, eM_slope, ncol=2)
ggsave('outputs//k600//roughness_theory//eM_channel.jpg', eM_channel, width=15, height=6)
