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

munge <- 0 #set to 1 if you need to generate results and individual river plots, otherwise 0

#constants-------------------
molarMass <- 44.01 #g/mol
t <- 25 #C
Sc_co2 <- 1742 + (-91.24*t) + (2.208*t^2) + (-0.0219*t^3) #Raymond etal 2012
pCO2_a <- 390 #uatm
henrys_law <- -60.2409+93.4517*(100/(t+273.15))+23.3585*log((t+273.15)/100) #Weiss 1974
g <- 9.8 #m/s2
sa <- 0 #m2 total surface area of 22 rivers

#legend labels------------------
legend_labels <- c('BIKER', 
                   'Model Using Observed velocity')

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
  BIKER_results <- read.csv('inputs/validation/results_SWOT_11day_markK.csv')
  BIKER_results$river <- substr(as.character(BIKER_results$river), 16, nchar(as.character(BIKER_results$river)))
  
  #Loop through SWOT rivers and just use Beaulieu etal 2012 CO2 data for all of them...------------------------
  output <- data.frame('river'=NA, 'timestep'=NA,
                       'FCO2_ulset'=NA, 
                       'FCO2_BIKER'=NA, 'FCO2_BIKER_low'=NA, 'FCO2_BIKER_high'=NA,
                       'FCO2_raymond2012'=NA, 'FCO2_raymond2012_low'=NA, 'FCO2_raymond2012_high'=NA,
                       'FCO2_raymond2013'=NA, 'FCO2_raymond2013_low'=NA, 'FCO2_raymond2013_high'=NA)
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
    reaches <- ncvar_get(data_in, 'River_Info/rch_bnd')
    
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
    
    #Calculate rating curves following raymond 2012 model
    V_raymond2012 <- exp(-1.64)*Q_obs^0.285
    V_raymond2012_low <- exp(-1.64+1.96*0.03)*Q_obs^(0.285+1.96*0.0091) #95% CIs
    V_raymond2012_high <- exp(-1.64-1.96*0.03)*Q_obs^(0.285-1.96*0.0091)
    
    #fit to the Raymond datapoints in Ulseth ds (marginally diff but I can get Cis)
    # V_raymond2012 <- exp(-1.84)*Q_obs^0.245
    # V_raymond2012_low <- exp(-1.9)*Q_obs^(0.218) #95% CIs
    # V_raymond2012_high <- exp(-1.78)*Q_obs^(0.272)

    #Calculate velocity using global rating curve from Liu etal in review
    V_raymond2013 <- exp(-1.06)*Q_obs^0.12
    V_raymond2013_low <- exp(-1.06)*Q_obs^(0.12+1.96*0.0091) #95% CIs but these are using 2012 uncertainities b/c there are no CIs published for these
    V_raymond2013_high <- exp(-1.06)*Q_obs^(0.12-1.96*0.0091)
    
    #eD
    eD <- g * V_obs * S_obs
    
    #Generate k600 estimates---------------------------
    k600_ulseth <- k600_craig(W_obs, S_obs, V_obs)
    k600_BIKER <- filter(BIKER_results, river == name)$kest_mean
    k600_raymond2012 <- k600_craig(W_obs, S_obs, V_raymond2012)
    k600_raymond2013 <- k600_craig(W_obs, S_obs, V_raymond2013)

    k600_BIKER_low <- filter(BIKER_results, river == name)$kest_low
    k600_BIKER_high <- filter(BIKER_results, river == name)$kest_high
    k600_raymond2012_low <- k600_craig(W_obs, S_obs, V_raymond2012_low)
    k600_raymond2012_high <- k600_craig(W_obs, S_obs, V_raymond2012_high)
    k600_raymond2013_low <- k600_craig(W_obs, S_obs, V_raymond2013_low)
    k600_raymond2013_high <- k600_craig(W_obs, S_obs, V_raymond2013_high)

    #Convert to kco2--------------------------
    kco2_ulseth <- colMeans((Sc_co2/600)^(-1/2)/k600_ulseth, na.rm=T)
    kco2_BIKER <- (Sc_co2/600)^(-1/2)/k600_BIKER
    kco2_raymond2012 <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012, na.rm=T)
    kco2_raymond2013 <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013, na.rm=T)

    kco2_BIKER_low <- (Sc_co2/600)^(-1/2)/k600_BIKER_low
    kco2_BIKER_high <- (Sc_co2/600)^(-1/2)/k600_BIKER_high
    kco2_raymond2012_low <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012_low, na.rm=T)
    kco2_raymond2012_high <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012_high, na.rm=T)
    kco2_raymond2013_low <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013_low, na.rm=T)
    kco2_raymond2013_high <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013_high, na.rm=T)

    #Calculate FCO2 using Beauliu data-------------------------------------------
    co2 <- Data$CO2_uatm[1:length(kco2_BIKER)] #measured co2
    
    FCO2_ulseth <- kco2_ulseth*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_BIKER <- kco2_BIKER*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2012 <- kco2_raymond2012*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2013 <- kco2_raymond2013*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

    FCO2_BIKER_low <- kco2_BIKER_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_BIKER_high <- kco2_BIKER_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2012_low <- kco2_raymond2012_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2012_high <- kco2_raymond2012_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2013_low <- kco2_raymond2013_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
    FCO2_raymond2013_high <- kco2_raymond2013_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

    #plot------------------------
    for_plot <- data.frame('river' = name, 'timestep'=1:length(FCO2_BIKER), 
                           'FCO2_ulset'=FCO2_ulseth, 
                           'FCO2_BIKER'=FCO2_BIKER, 'FCO2_BIKER_low'=FCO2_BIKER_low, 'FCO2_BIKER_high'=FCO2_BIKER_high,
                           'FCO2_raymond2012'=FCO2_raymond2012, 'FCO2_raymond2012_low'=FCO2_raymond2012_low, 'FCO2_raymond2012_high'=FCO2_raymond2012_high,
                           'FCO2_raymond2013'=FCO2_raymond2013, 'FCO2_raymond2013_low'=FCO2_raymond2013_low, 'FCO2_raymond2013_high'=FCO2_raymond2013_high)
    for_plot2 <- gather(for_plot, key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
     # gather(key=key_low, value=value_low, c(FCO2_ulset_low, FCO2_BIKER_low)) %>%
    #  gather(key=key_high, value=value_high, c(FCO2_ulset_high, FCO2_BIKER_high)) %>%
     # filter(substring(key, 6, 10) == substring(key_low, 6, 10) & substring(key, 6, 10) == substring(key_high, 6, 10))

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
    
    #ticker for total number of timesteps and total surface area
    sa <- sa + (mean(W_obs, na.rm = T) * reaches[length(reaches)]) #m2
    print(sa * 1e-6) #skm
  }
  
  output <- output[-1,]
  output$sa <- sa
  write.csv(output, 'outputs/flux_implications/CO2_results.csv')
}

#Summary stats-------------------------------
output <- read.csv('outputs/flux_implications/CO2_results.csv')
output$river <- as.character(output$river)
sa <- output[1,]$sa
output <- output[,-14]

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#total mass fluxes of CO2 from rivers------------------------------------------------
massFluxes <- gather(output, key=key, value=value, c(FCO2_ulset, FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013)) %>%
  group_by(key) %>%
  summarise(sumFCO2 = (sum(value, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365) #gG-C/yr from all rivers

masFlux_BIKER_low <- (sum(output$FCO2_BIKER_low, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365 #gG-C/yr from all rivers
masFlux_BIKER_high <- (sum(output$FCO2_BIKER_high, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365
masFlux_raymond2012_low <- (sum(output$FCO2_raymond2012_low, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365
masFlux_raymond2012_high <- (sum(output$FCO2_raymond2012_high, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365
masFlux_raymond2013_low <- (sum(output$FCO2_raymond2013_low, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365
masFlux_raymond2013_high <- (sum(output$FCO2_raymond2013_high, na.rm=T)/nrow(output)) * sa * 1e-9 * (12.011/44.01) * 365

massFluxes$low <- c(masFlux_BIKER_low, masFlux_raymond2012_low, NA, NA)
massFluxes$high <- c(masFlux_BIKER_high, masFlux_raymond2012_high, NA, NA)

#plot results----------------------------------------------------------------------
#Comparing velocity models----------------------------------------------------------------
#bar plots of mass fluxes
barPlots <- ggplot(massFluxes, aes(y=sumFCO2, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
 # geom_linerange(aes(ymin=low, ymax=high), size=1.2) +
  ylab('Bulk C Efflux [gG-C/yr]') +
  xlab('Velocity Model') +
  scale_fill_brewer(palette = 'Set1', name='', labels=c('BIKER', 'Raymond 2012', 'Raymond 2013', 'Observed')) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#CDFs comparing gauged rating curves versus our ungauged approach------------------------
t <- gather(output, key=key, value=value, c(FCO2_BIKER, FCO2_ulset, FCO2_raymond2012, FCO2_raymond2013))
t$flag <- ifelse(t$key == 'FCO2_ulset', 1, 0)
FCO2_cdfs <- ggplot(t, aes(x=value, color=key, linetype=factor(flag))) +
  stat_ecdf(size=1.25) +
  scale_color_brewer(palette='Set1') +
  xlab('FCO2 [g/m2*dy]') +
  ylab('Percentile') +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_hline(yintercept = 0.50, size=1.2) +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

FCO2_models <- plot_grid(barPlots, FCO2_cdfs, ncol=2, labels='auto')
ggsave('outputs/flux_implications/FCO2_models.jpg', FCO2_models, width=12, height=7)

#error stats and plot across all rivers and timesteps for BIKER vs obsereved V----------------------------------
lm <- lm(log(FCO2_BIKER)~log(FCO2_ulset), data=output)
lmr2 <- round(summary(lm)$r.squared,2)
rmse <- round(exp(Metrics::rmse(log(output$FCO2_ulset), log(output$FCO2_BIKER))), 2) #g/m2/dy

predInts <- predict(lm, interval='prediction')
output <- cbind(output, predInts)

flux_plot<- ggplot(output, aes(x=(FCO2_ulset), y=(FCO2_BIKER), color=key)) +
  geom_pointrange(aes(ymin = FCO2_BIKER_low, ymax = FCO2_BIKER_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='grey', method='lm', se=F)+
  geom_line(aes(y=exp(lwr)), color='grey', linetype='dashed', size=1.75) +
  geom_line(aes(y=exp(upr)), color='grey', linetype='dashed', size=1.75) +
  geom_abline(size=2, linetype='dashed', color='black') +
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
t <- filter(output, river=='Kanawha') %>%
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
Kanawha_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('Kanawha')

#Sacramento downstream
t <- filter(output, river=='OhioSection2') %>% 
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
OhioSection2_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('OhioSection2')

#Seine Upstream
t <- filter(output, river=='MississippiDownstream') %>% 
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
MississippiDownstream_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high))+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('MississippiDownstream')

#bring it together
legend <- get_legend(Kanawha_plot + theme(legend.box.margin = margin(0, 3, 3, 5)))
plot_fin <- plot_grid(Kanawha_plot + theme(legend.position = 'none'), OhioSection2_plot, MississippiDownstream_plot, legend, ncol=2, 
                      labels=c('b', 'c', 'd', NA))
plot_fin <- plot_grid(flux_plot, plot_fin, ncol=2, labels=c('a', NA))
ggsave('outputs/flux_implications/FCO2_plot.jpg', plot_fin, width=12, height=7)




#By river stats--------------------------------------------------------
for_plot <- gather(output, key=key, value=value, c(FCO2_raymond2012, FCO2_raymond2013, FCO2_BIKER))
stats_by_reach <- group_by(for_plot, river, key) %>%
  summarise(kge = hydroGOF::KGE(value, FCO2_ulset),
            nrmse = sqrt(mean((FCO2_ulset - value)^2)) / mean(FCO2_ulset, na.rm=T),
            rBIAS =   mean(value - FCO2_ulset) / mean(FCO2_ulset, na.rm=T),
            rrmse =   sqrt(mean((value - FCO2_ulset)^2 / FCO2_ulset^2)))

plot_stats <- gather(stats_by_reach, key=key2, value=value, c('nrmse', 'rBIAS', 'rrmse', 'kge'))
plotRivs <- ggplot(plot_stats, aes(x=key2, y=value, fill=key)) +
  geom_boxplot(size=1) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=-0.41, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  coord_cartesian(ylim = c(-1,1))+
  scale_fill_brewer(palette = 'Dark2', name='') +
  theme(legend.position = "bottom",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))
ggsave('outputs/flux_implications/FCO2_by_river.jpg', plotRivs, width=10, height=7)
