###########################
#Description: Modeling k600 from hydraulics/roughness. Using open-channel flow to explain k600 scaling dyanmics with surface turbulence
#Creator: Craig Brinkerhoff
#Spring 2021
###########################

print('starting k600 theory...')

options(scipen = 999)

#Raymond 2012 morphology stuff----------------------------------------------------
raymond2012 <- read.csv('inputs/k600/raymond2012.csv')
raymond2012 <- filter(raymond2012, is.na(k600)==0 & k600 > 0)
raymond2012$Reach.Control <- as.character(raymond2012$Reach.Control)
raymond2012$Bottom.Description <- as.character(raymond2012$Bottom.Description)

raymond2012$Reach.Control <- ifelse(raymond2012$Reach.Control != 'Channel Control' & raymond2012$Reach.Control != 'Pool and Riffle', 'Hybrid', raymond2012$Reach.Control)

raymond_summary <- group_by(raymond2012, Reach.Control) %>%
  tally()

qual_plot <- ggplot(raymond2012, aes(x=Reach.Control, fill=Reach.Control, y=k600)) +
  geom_violin(size=1.2, color='black', draw_quantiles=c(0.50)) +
  scale_fill_brewer(palette = 'Set2', name='') +
  geom_text(data = raymond_summary, aes(Reach.Control, Inf, label = n), vjust = 1, size=8)+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  xlab('') +
  ylab('K600 [m/day]') +
  theme(legend.position = 'none',
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs//k600//roughness_theory//qualPlot.jpg', qual_plot, width=9, height=7)

#Ulseth etal 2019------------------------------------------------------------
ulseth_data <- read.csv('data/k600/Ulseth_etal_2019.csv', fileEncoding="UTF-8-BOM") #contains data from 5 studies
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
data$keulegan <- 11*data$depth*(1/exp(data$Vms/(2.5*(9.8*data$depth*data$slope)^(1/2)))) #[m]

#log transform some variables-----------------------------------------------------------
data$log_eD <- log(data$eD)
data$log_k600 <- log(data$k600)
data$log_slope <- log(data$slope)
data$log_area <- log(data$area)

#partitoning of eD and eS---------------------------------------------------------
#Similar to Brinkerhoffetal 2019 function. Calibrating a Shield's grain size using Manning's equation, Shield's parameter, and Henderson relation f~(D/depth)^(1/3)
  #See math notes for the two velocity derivations and the math solving for the Henderson coefficient multiple ways. We ultimately use Cf=0.100
grainSizer <- function(v,s, d,cf, Sg){
  D <- 10^(seq(-5, 2, 0.005)) #virtual grain size [m] Tests 601 values

  #first method
  d_prime <- v^(3/2) * (sqrt((8*g*s)/(cf)))^(-3/2) * (1/D^(1/6))^(-3/2) #[m] #eq 15
  d_double_prime <- d - d_prime #[m]
  inverse_shields_prime <- ((Sg-1)*D)/(s*d_prime) #eq 16

  C_1 <- ifelse(inverse_shields_prime < 1.6, 41, 29)#eqs 20-22
                #ifelse(inverse_shields_prime < 17.8, 29, 23))
  C_2 <- ifelse(inverse_shields_prime < 1.6, 1.26,0.5)
                #ifelse(inverse_shields_prime < 17.8, 0.5, 0.415))
  v_m_double_prime_1 <- C_1/inverse_shields_prime^C_2 #eq 19, approximation of Einstein & Barbarossa sediment regimes

  #2nd method
  v_m_double_prime_2 <- v / sqrt(g*d_double_prime*s) #eq 14, closed form solution of v/m''

  #get calibrated value
  diff <- abs(v_m_double_prime_1-v_m_double_prime_2) #minimal difference b/w skin friction depth estimates
  D_fin <- D[which.min(diff)] #[m]

  #percent difference between estimates
  perc <- diff/v_m_double_prime_2
  perc <- min(perc, na.rm=T)*100

  return(c(D_fin, perc))
}

#run grain sizer
grainResults <- mapply(grainSizer,data$Vms,data$slope, data$depth, cf, Sg)
data$D35 <- grainResults[1,]
data$grain_perc_diff <- grainResults[2,]
summary(data$D35)
summary(data$grain_perc_diff)

data$log_D35 <- log(data$D35)

#calculate resistance partition coefficient using calibrated D35 (and its associated d')
data$inverse_shields_prime <- data$Vms^(3/2) * (cf/(8*g*data$D35))^(3/4) * data$slope^(1/4) * (Sg-1)^(-1)
data$d_prime <- data$inverse_shields_prime * data$slope^(-1) * (Sg-1) * data$D35 #calculated via shield's parameter and so D35. Remeber error is less than 5% so we good

data$f_prime <- cf *(data$D35/data$d_prime)^(1/3) #Henderson 1966

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

#3-piece eD model for BIKER----------------------------------------------------------------------
linearModel_eD <- lm(log_k600~log_eD, data=filter(data, eM_regime_craig == 'Medium'))
print(summary(linearModel_eD))
linearModel_eD <- lm(log_k600~log_eD, data=data)
lm_craig_eD <- segmented(linearModel_eD, npsi = 2) # two breakpoints
r2_craig_eD <- round(summary(lm_craig_eD)$r.squared, 2)
data$k600_eD <- exp(predict(linearModel_eD, data))

print(summary(lm_craig_eD))
print(slope(lm_craig_eD))
print(intercept(lm_craig_eD))

#2 piece D35 fit
linearModel_D35 <- lm(log_k600~log_D35, data=data)
lm_craig_D35 <- segmented(linearModel_D35, npsi = 1) # two breakpoints
r2_craig_D35 <- round(summary(lm_craig_D35)$r.squared, 2)
summary(lm_craig_D35)
data$k600_D35 <- exp(predict(lm_craig_D35, data))

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
ggsave('cache/k600/roughness_theory/eD_eS.jpg', eD_eS, width=9, height=7)

data$depthFlag <- ifelse(data$depth < 10, 'Roughness-based', 'Wind-based')
eM_k600 <- ggplot(data, aes(x=eM, y=k600, color=depthFlag)) +
  geom_point(size=4, alpha=0.50) +
  geom_line(aes(x=eM, y=k600_eM), color='darkgreen', size=3) +
  geom_hline(yintercept=35, color='black', linetype='dashed', size=2)+
  annotate("text", label = "Ulseth: ~maximum \ndiffusive gas exchange (35 m/dy)", x = 10^-3, y = 10^2, size = 5, colour = "black") +
  xlab('eM [W/kg]') +
  ylab('k600 [m/dy]') +
  annotate("text", label = paste0('r2: ', r2_craig_eM), x = 10^-4, y = 10^3, size = 5, colour = "darkgreen")+
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
ggsave('cache/k600/roughness_theory/eM_k600.jpg', eM_k600, width=9, height=7)

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
ggsave('cache/k600/roughness_theory/eM_partition.jpg', eM_partition, width=9, height=7)

#ulseth_line <- data.frame('D'=seq(0.02, 0.1, 0.01))
#ulseth_line$k600 <- exp(3.41+30.52*ulseth_line$D)

eM_D35 <- ggplot(data, aes(x=D35, y=k600, color=slope)) +
  geom_point(size=4, alpha=0.75) +
  geom_line(aes(x=D35, y=k600_D35), color='darkgreen', size=3) +
#  geom_line(data=ulseth_line, aes(x=D, y=k600), color='darkblue', size=2)+ #Ulseth roughness model
  geom_hline(yintercept=35, color='black', linetype='dashed', size=2)+
  annotate("text", label = "Ulseth: ~maximum \ndiffusive gas exchange (35 m/dy)", x = 10^-3, y = 10^2, size = 5, colour = "black") +
  xlab('Bed D35 [m]') +
  ylab('k600 [m/dy]') +
  annotate("text", label = paste0('r2: ', r2_craig_D35), x = 10^-4, y = 10^3, size = 5, colour = "darkgreen")+
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
ggsave('cache/k600/roughness_theory/eM_D35.jpg', eM_D35, width=9, height=7)

partitionPlot <- ggplot(data, aes(x=partition*100, color=eM_regime_craig)) +
    stat_ecdf(size=2) +
    xlab('Turbulence due to \nDepth-Scale Form Drag [%]') +
    scale_color_brewer(palette='Set2', name='eM breakpoint regime')+
    theme(axis.text=element_text(size=20),
      axis.title=element_text(size=24,face="bold"),
      legend.text = element_text(size=17),
      legend.title = element_text(size=17, face='bold'),
     legend.position='none')
ggsave('cache/k600/roughness_theory/partitionPlot.jpg', partitionPlot, width=9, height=7)

eM_area <- ggplot(data, aes(x=eM, y=area)) +
  geom_point(size=4, alpha=0.50) +
  geom_line(aes(x=eM, y=area_pred), color='darkgreen', size=3) +
  xlab('eM [W/kg]') +
  ylab('Channel Area [m2]') +
  annotate("text", label = paste0('r2: ', r2_craig_area), x = 10^-4, y = 10^0, size = 5, colour = "darkgreen")+
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
  annotate("text", label = paste0('r2: ', r2_craig_slope), x = 10^-4, y = 10^0, size = 5, colour = "darkgreen")+
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
ggsave('cache/k600/roughness_theory/eM_channel.jpg', eM_channel, width=15, height=6)


#k600 prior for BIKER--------------------------
#Hagemann et al. (2017)-style prior using slope. IMPLEMENTED.
#lmSlope <- lm(log_k600~log_slope, data=filter(data, eM_regime_craig == 'Medium'))
#print(summary(lmSlope))
#segLmSlope <- segmented(lmSlope, npsi=2) #r2 0.55
#slope(segLmSlope)
#intercept(segLmSlope)
#summary(segLmSlope)

slopePlot <- ggplot(filter(data, eM_regime_craig == 'Medium'), aes(x=log_slope, y=log_k600)) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('cache/k600/roughness_theory/slopeK600.jpg', slopePlot, width=15, height=6)
