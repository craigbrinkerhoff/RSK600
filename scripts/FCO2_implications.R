#Toy model to explore implications of k model for GHG fluxes
  #various upscaling models versus BIKER algorithm
  #same script as FCH4_implications but for carbon dioxide
#Craig Brinkerhoff
#Winter 2020

library(tidyverse)
library(ncdf4)
library(cowplot)
library(ggtext)
theme_set(theme_cowplot())

munge <- 1 #set to 1 if you need to generate results and individual river plots, otherwise 0

#constants-------------------
molarMass <- 44.01 #g/mol
t <- 25 #C
Sc_co2 <- 1742 + (-91.24*t) + (2.208*t^2) + (-0.0219*t^3) #Raymond etal 2012
pCO2_a <- 390 #uatm
henrys_law <- -60.2409+93.4517*(100/(t+273.15))+23.3585*log((t+273.15)/100) #Weiss 1974
g <- 9.8 #m/s2

#legend labels------------------
legend_labels <- c('BIKER', 
                   # '\nRaymond Model Used Globally in \n1) Raymond etal 2013 \n2) Lauerwald etal 2015 \n3) Liu etal in review \n',
                   # 'This study',
                   'RS-able Ulseth etal 2019')

#Functions----------------------------------------
calcdA_mat <- function(w, h) {
  stopifnot(all(dim(w) == dim(h)))
  dA <- w
  for (i in 1:nrow(dA)) {
    dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
  }
  
  dA
}

calcdA_vec <- function(w, h) {
  words <- order(w)
  warr <- w[words]
  harr <- h[words]
  delh <- c(0, diff(harr))
  delA <- cumsum(warr * delh)
  dA <- 1:length(w)
  dA[words] <- delA
  dA
}

k600_craig <- function(w, s, v) {
  widthRegime <- ifelse(mean(w, na.rm = T) < 10 & mean(s) < 0.05, 1,
                        ifelse(mean(w, na.rm = T) < 10 & mean(s, na.rm=T) >= 0.05, 5,
                               ifelse(mean(w, na.rm = T) < 50, 2,
                                      ifelse(mean(w, na.rm = T)< 100, 3, 4))))
  if (widthRegime == 1){
    k <- (111.58121 * (g*s * (v))^0.6488131)
  }
  if(widthRegime ==2 ){
    k <- (109.04977 * (g*s * (v))^0.6624354)
  }
  if(widthRegime == 3){
    k <- (31.84344 * (g*s * (v))^0.4366114)
  }
  if(widthRegime == 4) {
    k <- (14.16939 * (g*s * (v))^0.2724834)
  }
  if(widthRegime == 5) {
    k <- (792.63149 * (g*s * (v))^1.3160065)
  }
  return(k)
}

#run-------------------------------------------
if (munge == 1) {
  #read in co2 data (Beauliu etal 2012)--------------------------------
  data <- read.csv('inputs/flux_implications/Beauliu2012.csv', header = TRUE)
  data <- data[2:30,]
  Data <- select(data, c(Date, CO2))
  #Data$Date <- as.Date(as.character(Data$Date), '%m/%d/%Y')
  colnames(Data) <- c('Date', 'CO2_umol_L')
  Data$CO2_umol_L <- as.numeric(as.character(Data$CO2_umol_L))
  Data$CO2_uatm <- Data$CO2_umol_L / exp(henrys_law)
  
  #read in BIKER results-----------------------------------------------
  BIKER_results <- read.csv('inputs/validation/results_SWOT_11day_slopeK.csv')
  BIKER_results$river <- substr(as.character(BIKER_results$river), 16, nchar(as.character(BIKER_results$river)))
  
  #Loop through SWOT rivers and just use Beaulieu etal 2012 CO2 data for all of them...------------------------
  output <- data.frame('river'=NA, 'timestep'=NA,
                       'FCO2_raymo'=NA, 'FCO2_raymo_low'=NA, 'FCO2_raymo_high'=NA, 
                       'FCO2_ulset'=NA, 'FCO2_ulset_low'=NA, 'FCO2_ulset_high'=NA, 
                       'FCO2_BIKER'=NA, 'FCO2_BIKER_low'=NA, 'FCO2_BIKER_high'=NA)
  files = list.files('inputs/flux_implications', pattern="*.nc", full.names = TRUE)
  for (river in files[2:length(files)]){ #skip Arial Khan
    #read in swot rivers
    name <- substr(river, 26, nchar(river))
    name <- substr(name,1,nchar(name)-3)
    data_in =nc_open(river)
    
    W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
    H_obs=ncvar_get(data_in,'Reach_Timeseries/H')
    S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
    area = ncvar_get(data_in, 'Reach_Timeseries/A')
    Q_obs = ncvar_get(data_in, 'Reach_Timeseries/Q')
    
    #prep pepsi rivers
    S_obs[S_obs<=0]=NA
    S_obs[is.na(S_obs)] = 0.000001 #min obs SWOT slope
    W_obs[W_obs<0]=NA
    H_obs[H_obs<0]=NA
    area[area<0]=NA
    Q_obs[Q_obs<0]=NA
    
    #Sample every 11 days (swot-style)
    W_obs = W_obs[,seq(1, ncol(W_obs), 11)]
    S_obs = S_obs[,seq(1, ncol(S_obs), 11)]  
    H_obs = H_obs[,seq(1, ncol(H_obs), 11)]  
    Q_obs = Q_obs[,seq(1, ncol(Q_obs), 11)]  
    area = area[,seq(1, ncol(area), 11)]
    
    if (any(apply(S_obs,2,sum,na.rm=TRUE) ==0)){
      remove_index =  which((apply(S_obs,2,sum,na.rm=TRUE) ==0) ==TRUE)
      
      W_obs=W_obs[,-remove_index]
      H_obs=H_obs[,-remove_index]
      S_obs=S_obs[,-remove_index]
      Q_obs=Q_obs[,-remove_index]
      area=area[,-remove_index]
    }
    
    if (any(apply(S_obs,1,sum,na.rm=TRUE) ==0)){
      remove_index =  which((apply(S_obs,1,sum,na.rm=TRUE) ==0) ==TRUE)
      
      W_obs=W_obs[-remove_index,]
      H_obs=H_obs[-remove_index,]
      S_obs=S_obs[-remove_index,]
      Q_obs=Q_obs[,-remove_index]
      area=area[,-remove_index]
    }
    
    #Calculate dA matrix from RS W and H
    dA_obs <- calcdA_mat(W_obs,H_obs) #m2
    
    #calculate observed average flow velocity
    V_obs <- Q_obs / area #m/s
    
    #Calculate k600 and kco2 across models------------------------------------
    eD <- V_obs * S_obs * g #m2/s3
    
    #Generate k600 estimates---------------------------
    k600_raymond <- V_obs*S_obs*2841 + 2.02 #Raymond etal 2012 eq 4 (used in Raymond etal 2013 / Lauerwald etal 2015 / Liu etal in review)
   # k600_ulseth <- exp(ifelse(eD < 0.02, 3.10+0.35*log(eD), 6.43+1.18*log(eD))) #ulseth etal 2019
    k600_ulseth <- k600_craig(W_obs, S_obs, V_obs)
    k600_BIKER <- filter(BIKER_results, river == name)$kest_mean
    
    k600_BIKER_low <- filter(BIKER_results, river == name)$kest_low
    k600_BIKER_high <- filter(BIKER_results, river == name)$kest_high
    k600_raymond_low <- V_obs*S_obs*(2841-107) + (2.02-0.209) 
    k600_raymond_high <- V_obs*S_obs*(2841+107) + (2.02+0.209) 
    k600_ulseth_low <- exp(ifelse(eD < 0.02, 3.10+0.31*log(eD), 6.43+1.10*log(eD)))
    k600_ulseth_high <- exp(ifelse(eD < 0.02, 3.10+0.41*log(eD), 6.43+1.30*log(eD)))
    
    #Convert to kco2--------------------------
    kco2_raymond <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond, na.rm = T)
    kco2_ulseth <- colMeans((Sc_co2/600)^(-1/2)/k600_ulseth, na.rm=T)
    kco2_BIKER <- (Sc_co2/600)^(-1/2)/k600_BIKER

    kco2_BIKER_low <- (Sc_co2/600)^(-1/2)/k600_BIKER_low
    kco2_BIKER_high <- (Sc_co2/600)^(-1/2)/k600_BIKER_high
    kco2_raymond_low <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond_low, na.rm = T)
    kco2_raymond_high <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond_high, na.rm = T)
    kco2_ulseth_low <- colMeans((Sc_co2/600)^(-1/2)/k600_ulseth_low, na.rm = T)
    kco2_ulseth_high <- colMeans((Sc_co2/600)^(-1/2)/k600_ulseth_high, na.rm = T)
    
    #Calculate FCO2 using Beauliu data-------------------------------------------
    co2 <- Data$CO2_uatm[1:length(kco2_BIKER)] #measured co2
    
    FCO2_raymond <- kco2_raymond*((co2-pCO2_a)*exp(henrys_law)*molarMass*(1/1000000)*(1/0.001)) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_ulseth <- kco2_ulseth*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_BIKER <- kco2_BIKER*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

    FCO2_BIKER_low <- kco2_BIKER_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_BIKER_high <- kco2_BIKER_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond_low <- kco2_raymond_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond_high <- kco2_raymond_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_ulseth_low <- kco2_ulseth_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_ulseth_high <- kco2_ulseth_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    
    #plot------------------------
    for_plot <- data.frame('river' = name, 'timestep'=1:length(FCO2_BIKER), 
                           'FCO2_raymo'=FCO2_raymond, 'FCO2_raymo_low'=FCO2_raymond_low, 'FCO2_raymo_high'=FCO2_raymond_high, 
                           'FCO2_ulset'=FCO2_ulseth, 'FCO2_ulset_low'=FCO2_ulseth_low, 'FCO2_ulset_high'=FCO2_ulseth_high, 
                           'FCO2_BIKER'=FCO2_BIKER, 'FCO2_BIKER_low'=FCO2_BIKER_low, 'FCO2_BIKER_high'=FCO2_BIKER_high)
    for_plot2 <- gather(for_plot, key=key, value=value, c(FCO2_ulset, FCO2_BIKER)) %>%
      gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
      gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
      filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10))
    
    plot <- ggplot(for_plot2, aes(x=timestep, y=value, color=key)) +
      geom_line(size=1.5) +
     # geom_pointrange(aes(ymin=value_low, ymax=value_high))+
      ylab('FCO2 [g/m2*dy]') +
      xlab('Timestep') +
      theme(legend.position = 'bottom') +
      ggtitle(name)+
      scale_color_brewer(palette='Set2', name='', labels=legend_labels)
    ggsave(paste0('outputs/flux_implications/by_river/FCO2_', name, '.jpg'), plot, width=7, height=7)
  
    #write results to file----------------
    output <- rbind(output, for_plot)
  }
  output <- output[-1,]
  write.csv(output, 'outputs/flux_implications/CO2_results.csv')
}

#Summary stats-------------------------------
output <- read.csv('outputs/flux_implications/CO2_results.csv')
output$river <- as.character(output$river)

#taking care of 1 thing:
  #1) CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#plot results-----------------------------------------------------
output2 <- gather(output, key=key, value=value, c(FCO2_ulset, FCO2_BIKER)) %>%
  gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
  gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
  filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10)) %>%
  group_by(timestep, key) %>%
  summarise(sumFCO2 = sum(value, na.rm=T),
            sumFCO2_low = sum(value_low, na.rm=T),
            sumFCO2_high = sum(value_high, na.rm=T))

#error stats across all rivers and timesteps
lm <- lm(log(FCO2_BIKER)~log(FCO2_ulset), data=output)
lmr2 <- round(summary(lm)$r.squared,2)
rmse <- round(exp(Metrics::rmse(log(output$FCO2_ulset), log(output$FCO2_BIKER))), 2) #g/m2/dy

predInts <- predict(lm, interval='prediction')
output <- cbind(output, predInts)

flux_plot<- ggplot(output, aes(x=(FCO2_ulset), y=(FCO2_BIKER))) +
  geom_point(size=5, pch=21, color='black', fill='#1b9e77') +
  geom_smooth(size=2, color='black', method='lm', se=F)+
  geom_line(aes(y=exp(lwr)), color='black', linetype='dashed', size=1.75) +
  geom_line(aes(y=exp(upr)), color='black', linetype='dashed', size=1.75) +
  geom_abline(size=2, linetype='dashed', color='grey') +
  xlab('FCO2 via observed velocity [g/m2*dy]') +
  ylab('FCO2 via BIKER [g/m2*dy]') +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold')) +
  geom_richtext(aes(x=10^-0.4, y=10^0.5), color='black', label=paste0('RMSE: ', rmse, ' g/m2*dy')) +
  geom_richtext(aes(x=10^-0.5, y=10^0.4), color='black', label=paste0('r<sup>2</sup>: ', lmr2), size=5)

#add some individual river plots
#Cumberland
t <- filter(output, river=='Cumberland') %>%
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER)) %>%
  gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
  gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
  filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10))
Cumberland_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  #geom_pointrange(aes(ymin=value_low, ymax=value_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels)

#Sacramento downstream
t <- filter(output, river=='SacramentoDownstream') %>% 
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER)) %>%
  gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
  gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
  filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10))
SacDown_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
 # geom_pointrange(aes(ymin=value_low, ymax=value_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels)

#Seine Upstream
t <- filter(output, river=='OhioSection5') %>% 
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER)) %>%
  gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
  gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
  filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10))
Ohio5_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
#  geom_pointrange(aes(ymin=value_low, ymax=value_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels)

#bring it together
legend <- get_legend(Cumberland_plot + theme(legend.box.margin = margin(0, 3, 3, 5)))
plot_fin <- plot_grid(Cumberland_plot + theme(legend.position = 'none'), SacDown_plot, Ohio5_plot, legend, ncol=2, 
                      labels=c('Cumberland', 'Sacramento\nDownstream', 'Ohio\nSection 5', NA))
plot_fin <- plot_grid(flux_plot, plot_fin, ncol=2)
ggsave('outputs/flux_implications/FCO2_plot.jpg', plot_fin, width=12, height=7)




#scraps-------------------------------------------------
# flux_plot <- ggplot(output2, aes(x=timestep, y=sumFCO2, color=key)) +
#   geom_line(size=1.25, linetype='dashed') +
#   geom_point(size=5)+
# #  geom_pointrange(aes(ymin=sumFCO2_low, ymax=sumFCO2_high)) +
#   ylab('FCO2 [g/m2*dy]') +
#   xlab('Timestep') +
#   annotate(x=25.5, y=60, 'text', label='Timesteps \ndo not necessairly \ninclude all rivers \nbecause of variable \ntimeseries lengths', size=4)+
#   annotate(x=25.5, y=45, 'text', label='k600 = a_i(gSV)^b_i', size=4, color='#fc8d62')+ #e78ac3
#  # annotate(x=25.5, y=42.5, 'text', label='k600 = a(gSV)^b', size=4, color='#8da0cb')+
#   annotate(x=25.5, y=39, 'text', label='Bayesian inference\nof a_i(gSV)^b_i', size=4, color='#66c2a5')+
#  # annotate(x=25.5, y=35, 'text', label='k600 = B0 + B1(SV)', size=4, color='#fc8d62')+
#   theme(legend.position = 'none') +
#   geom_richtext(aes(x=5, y=50), color='black', label=paste0('r<sup>2</sup>: ', lmr2), size=5) +
#   geom_richtext(aes(x=5, y=52), color='black', label=paste0('RMSE: ', rmse, ' g/m2*dy'), size=5) +
#   ggtitle('Upscaled Across 22 SWOT rivers')+
#   scale_color_brewer(palette='Set2', name='', labels=legend_labels)
# flux_plot