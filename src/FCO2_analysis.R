###########################
#Description: Uses BIKER and rating curves to calculate FCO2 and bulk C fluxes from the SWOT simulated rivers
#Craig Brinkerhoff
#Winter 2021
##########################

print('starting FCO2...')

#legend labels------------------
legend_labels <- c('BIKER',
                   'Model Using Observed velocity')

#Function to generate FCO2 from all velocity models and all swot rivers
fun_FCO2_analysis <- function(river) {
  #read in swot rivers
  name <- substr(river, 24, nchar(river))
  name <- substr(name,1,nchar(name)-3)
  data_in =nc_open(river)

  W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
  H_obs=ncvar_get(data_in,'Reach_Timeseries/H')
  S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
  area = ncvar_get(data_in, 'Reach_Timeseries/A')
  Q_obs = ncvar_get(data_in, 'Reach_Timeseries/Q')
  reaches <- ncvar_get(data_in, 'River_Info/rch_bnd')

  #Account for models known to not start on 1/1/XXXX. Seine starts in june and ___ starts in Feburary
  start_date <- ncvar_get(data_in, 'Reach_Timeseries/t')
  start_date <- ifelse(start_date[1] !=1, start_date[1], 1)
  river_start_date <- format(as.Date('0000/01/01') + start_date, format='%m')

  #prep pepsi rivers
  S_obs[S_obs<=0]=NA
  S_obs[is.na(S_obs)] = 0.000001 #min obs SWOT slope
  W_obs[W_obs<0]=NA
  H_obs[H_obs<0]=NA
  area[area<0]=NA
  Q_obs[Q_obs<0]=NA

  #Sample every 11 days (swot-style)
  W_obs = W_obs[,seq(river_start_date, ncol(W_obs), 11)]
  S_obs = S_obs[,seq(river_start_date, ncol(S_obs), 11)]
  H_obs = H_obs[,seq(river_start_date, ncol(H_obs), 11)]
  Q_obs = Q_obs[,seq(river_start_date, ncol(Q_obs), 11)]
  area = area[,seq(river_start_date, ncol(area), 11)]

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

  s <- rbind(s, data.frame('river'=river, meanS=mean(S_obs, na.rm=T)))

  #Calculate dA matrix from RS W and H
  dA_obs <- calcdA_mat(W_obs,H_obs) #m2

  #calculate observed average flow velocity
  V_obs <- Q_obs / area #m/s
  D_obs <- area / W_obs #m

  #Calculate rating curves following raymond 2012 model
  V_raymond2012 <- exp(-1.64)*Q_obs^0.285
  V_raymond2012_low <- exp(-1.64+1.96*0.03)*Q_obs^(0.285+1.96*0.0091) #95% CIs
  V_raymond2012_high <- exp(-1.64-1.96*0.03)*Q_obs^(0.285-1.96*0.0091)

  D_raymond2012 <- exp(-0.895)*Q_obs^0.294

  #Calculate velocity using global rating curve from 1/2 the raymond 2013 model (used in Lauerwald 2015)
  V_raymond2013 <- exp(-1.06)*Q_obs^0.12
  V_raymond2013_low <- exp(-1.06)*Q_obs^(0.12+1.96*0.0091) #95% CIs but these are using 2012 uncertainities b/c there are no CIs published for these
  V_raymond2013_high <- exp(-1.06)*Q_obs^(0.12-1.96*0.0091)

  D_raymond2013 <- exp(1/(1.86*-1.06))*Q_obs^(1-0.51-0.12)

  V_lauerwald2015 <- V_raymond2013
  V_lauerwald2015_low <- V_raymond2013_low
  V_lauerwald2015_high <- V_raymond2013_high

  D_lauerwald2015 <- D_raymond2013

  #avgerage of the models was the actual implementation in Raymond 2013
  d <- data.frame(V_raymond2012, V_raymond2013)
  V_raymond2013 <- rowMeans(d)

  d <- data.frame(D_raymond2012, D_raymond2013)
  D_raymond2013 <- rowMeans(d)

  #eD
  eD <- g * V_obs * S_obs

  #ustar
  Ustar <- sqrt(g*S_obs*D_obs)

  #Generate k600 estimates---------------------------
  k600_ulseth <- k600_craig(S_obs, V_obs)
  k600_BIKER <- filter(BIKER_results, river == name)$kest_mean
  k600_raymond2012 <- k600_craig(S_obs, V_raymond2012)
  k600_raymond2013 <- k600_craig(S_obs, V_raymond2013)
  k600_Lauerwald2015 <- k600_craig(S_obs, V_lauerwald2015)

  k600_BIKER_low <- filter(BIKER_results, river == name)$kest_low
  k600_BIKER_high <- filter(BIKER_results, river == name)$kest_high
  k600_raymond2012_low <- k600_craig(S_obs, V_raymond2012_low)
  k600_raymond2012_high <- k600_craig(S_obs, V_raymond2012_high)
  k600_raymond2013_low <- k600_craig(S_obs, V_raymond2013_low)
  k600_raymond2013_high <- k600_craig(S_obs, V_raymond2013_high)
  k600_Lauerwald2015_low <- k600_craig(S_obs, V_lauerwald2015_low)
  k600_Lauerwald2015_high <- k600_craig(S_obs, V_lauerwald2015_high)

  #Generate kL20 estimates---------------------------
  kL20_ulseth <- ko2_craig(D_obs, S_obs)
  kL20_BIKER <- filter(BIKER_results, river == name)$kest_mean
  kL20_BIKER <- kL20_BIKER[seq(1, length(kL20_BIKER), 11)]
  kL20_BIKER_low <- filter(BIKER_results, river == name)$kest_low
  kL20_BIKER_low <- kL20_BIKER_low[seq(1, length(kL20_BIKER_low), 11)]
  kL20_BIKER_high <- filter(BIKER_results, river == name)$kest_high
  kL20_BIKER_high <- kL20_BIKER_high[seq(1, length(kL20_BIKER_high), 11)]
  kL20_raymond2012 <- ko2_craig(D_raymond2012, S_obs)
  kL20_raymond2013 <- ko2_craig(D_raymond2013, S_obs)
  kL20_Lauerwald2015 <- ko2_craig(D_lauerwald2015, S_obs)

  #Calculate FCO2 and Sc using Beauliu data-------------------------------------------
  #sort by starting mmonth here
  Data_sort <- Data[(as.numeric(river_start_date)*2):nrow(Data),]
  co2 <- Data_sort$CO2_uatm[1:length(kL20_BIKER)] #measured co2
  Sc_co2 <- Sc_co2_func(Data_sort$Water_temp_C) #Sc from measured water temperature
  Sc_co2 <- Sc_co2[1:length(kL20_BIKER)]
  henrys_law <- henrys_law_func(Data_sort$Water_temp_C) #Henry's solubility coefficient from measured water temperature
  henrys_law <- henrys_law[1:length(kL20_BIKER)]

  #Convert to kco2--------------------------
#  kco2_ulseth <- colMeans((Sc_co2/600)^(-1/2)/k600_ulseth, na.rm=T)
#  kco2_BIKER <- (Sc_co2/600)^(-1/2)/k600_BIKER
#  kco2_raymond2012 <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012, na.rm=T)
#  kco2_raymond2013 <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013, na.rm=T)
#  kco2_Lauerwald2015 <- colMeans((Sc_co2/600)^(-1/2)/k600_Lauerwald2015, na.rm=T)

  kco2_BIKER_low <- (Sc_co2/600)^(-1/2)/k600_BIKER_low
  kco2_BIKER_high <- (Sc_co2/600)^(-1/2)/k600_BIKER_high
  kco2_raymond2012_low <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012_low, na.rm=T)
  kco2_raymond2012_high <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2012_high, na.rm=T)
  kco2_raymond2013_low <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013_low, na.rm=T)
  kco2_raymond2013_high <- colMeans((Sc_co2/600)^(-1/2)/k600_raymond2013_high, na.rm=T)
  kco2_Lauerwald2015_low <- colMeans((Sc_co2/600)^(-1/2)/k600_Lauerwald2015_low, na.rm=T)
  kco2_Lauerwald2015_high <- colMeans((Sc_co2/600)^(-1/2)/k600_Lauerwald2015_high, na.rm=T)

  kco2_ulseth <- colMeans((Sc_co2/530)^(-1/2)*kL20_ulseth, na.rm=T)
  kco2_BIKER <- (Sc_co2/530)^(-1/2)*kL20_BIKER
  kco2_raymond2012 <- colMeans((Sc_co2/530)^(-1/2)*kL20_raymond2012, na.rm=T)
  kco2_raymond2013 <- colMeans((Sc_co2/530)^(-1/2)*kL20_raymond2013, na.rm=T)
  kco2_Lauerwald2015 <- colMeans((Sc_co2/530)^(-1/2)*kL20_Lauerwald2015, na.rm=T)
  kco2_BIKER_low <- (Sc_co2/530)^(-1/2)*kL20_BIKER_low
  kco2_BIKER_high <- (Sc_co2/530)^(-1/2)*kL20_BIKER_high

  #Obtain actual CO2 fluxes--------------------------------------------
  FCO2_ulseth <- kco2_ulseth*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_BIKER <- kco2_BIKER*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2012 <- kco2_raymond2012*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2013 <- kco2_raymond2013*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_Lauerwald2015 <- kco2_Lauerwald2015*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

  FCO2_BIKER_low <- kco2_BIKER_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_BIKER_high <- kco2_BIKER_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2012_low <- kco2_raymond2012_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2012_high <- kco2_raymond2012_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2013_low <- kco2_raymond2013_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2013_high <- kco2_raymond2013_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_Lauerwald2015_low <- kco2_Lauerwald2015_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_Lauerwald2015_high <- kco2_Lauerwald2015_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

  #Surface area
  #first, get reach lengths from the cumulative reach lengths in river data
  r <- reaches
  for (i in 2:length(reaches)) {
    r[i] <- reaches[i] - reaches[i-1]
  }
  r <- r[-1]

  #then calculate surface area
  sa <- sum(mean(W_obs, na.rm = T) * r) #m2

  #plot------------------------
  for_plot <- data.frame('river' = name, 'timestep'=1:length(FCO2_BIKER), 'sa'=sa,
                        'FCO2_ulset'=FCO2_ulseth,
                        'FCO2_BIKER'=FCO2_BIKER, 'FCO2_BIKER_low'=FCO2_BIKER_low, 'FCO2_BIKER_high'=FCO2_BIKER_high,
                        'FCO2_raymond2012'=FCO2_raymond2012, 'FCO2_raymond2012_low'=FCO2_raymond2012_low, 'FCO2_raymond2012_high'=FCO2_raymond2012_high,
                        'FCO2_raymond2013'=FCO2_raymond2013, 'FCO2_raymond2013_low'=FCO2_raymond2013_low, 'FCO2_raymond2013_high'=FCO2_raymond2013_high,
                        'FCO2_Lauerwald2015'=FCO2_Lauerwald2015, 'FCO2_Lauerwald2015_low'=FCO2_Lauerwald2015_low, 'FCO2_Lauerwald2015_high'=FCO2_Lauerwald2015_high)

    #write results to file----------------
    write.csv(for_plot, paste0('cache/FCO2/by_river/results_', name, '.csv'))
}

s <- data.frame('river'=NA, 'meanS'=NA)
#run-------------------------------------------
#read in co2 and water temperature data (Beauliu etal 2012)--------------------------------
data <- read.csv('data/FCO2/Beauliu2012.csv', header = TRUE)
data <- data[5:30,] #only keep 26 biweekly sample so they reflect a single year (fall 2009-winter/spring/summer 2008). This amounts to
Data <- select(data, c(Date, CO2, Temperature))
Data$Date <- as.Date(as.character(Data$Date), '%m/%d/%y')
Data$Date <- format(Data$Date, format="%m-%d")
Data <- Data[order(Data$Date),]
colnames(Data) <- c('Date', 'CO2_umol_L', 'Water_temp_C')
Data$Water_temp_C <- as.numeric(as.character(Data$Water_temp_C))
Data$CO2_umol_L <- as.numeric(as.character(Data$CO2_umol_L))
Data$CO2_uatm <- Data$CO2_umol_L / exp(henrys_law_func(Data$Water_temp_C))

#read in BIKER results (only those with no measurement error)-----------------------------------------------
BIKER_results <- read.csv('cache/validation/BIKER_validation_results.csv')
BIKER_results <- filter(BIKER_results, errFlag == 0)

#SWOT rivers (Frasson etal 2019 & Rodriquez etal 2020)----------------------------------------------------
files <- list.files('data/Frasson_etal_2019', pattern="*.nc", full.names = TRUE)
files <- files[-1] #remove Arial Khan

#Run parallelized function------------------------------
system.time(
  results <- mclapply(files, fun_FCO2_analysis, mc.cores=cores)
)

#results <- fun_FCO2_analysis(files[1])

#Make final full results file----------------------------
df_fin <- list.files(path='cache/FCO2/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows
write.csv(df_fin, 'cache/FCO2/CO2_results.csv')

#plot Beaulieu 2012 data--------------------------------
beaulieu <- ggplot(data=Data, aes(y=CO2_uatm, x=Date)) +
  geom_point(size=5, color='darkgreen') +
  ylab('Water-Side CO2 [uatm]')
ggsave('cache/FCO2/Beaulieu_timeseries.jpg', beaulieu, width=15, height=6)

#Summary stats-------------------------------
output <- read.csv('cache/FCO2/CO2_results.csv')
output$river <- as.character(output$river)

total_sa <- output[!duplicated(output$sa), ] #m2
total_sa <- sum(total_sa$sa)

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#total mass fluxes of CO2 from rivers (for barplots)------------------------------------------------
massFluxes <- gather(output, key=key, value=value, c(FCO2_ulset, FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013, FCO2_Lauerwald2015)) %>%
  group_by(key) %>%
  summarise(sumFCO2 = (sum(value, na.rm=T)/nrow(output)) * total_sa * 1e-9 * (12.011/44.01) * 365) #gG-C/yr from all rivers
write.csv(massFluxes, 'cache/FCO2/bulkFluxes.csv')

#plot results----------------------------------------------------------------------
#bar plots of mass fluxes
barPlots <- ggplot(massFluxes, aes(y=sumFCO2, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_ulset', ]$sumFCO2, size=1.2, linetype='dashed') +
  ylab('Bulk C Efflux [gG-C/yr]') +
  xlab('Velocity Model') +
  scale_fill_brewer(palette = 'Set1', name='', labels=c('BIKER', 'Lauerwald \n2015', 'Raymond \n2012', 'Raymond \n2013', 'Observed')) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#By river performance stats
for_plot <- gather(output, key=key, value=value, c(FCO2_raymond2012, FCO2_raymond2013, FCO2_BIKER, FCO2_Lauerwald2015))
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
  coord_cartesian(ylim = c(-1,1.5))+
  scale_fill_brewer(palette = 'Set1', name='') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

FCO2_models <- plot_grid(barPlots, plotRivs, ncol=2, labels='auto')
ggsave('cache/FCO2/FCO2_models.jpg', FCO2_models, width=12, height=6)

#error stats and plot across all rivers and timesteps
lm <- lm(log10(FCO2_BIKER)~log10(FCO2_ulset), data=output)
lmr2 <- round(summary(lm((FCO2_BIKER)~(FCO2_ulset), data=output))$r.squared,2)
rmse <- round((Metrics::rmse((output$FCO2_ulset), (output$FCO2_BIKER))), 2) #g/m2/dy

predInts <- predict(lm, interval='prediction')
output <- cbind(output, predInts)

flux_plot<- ggplot(output, aes(x=(FCO2_ulset), y=(FCO2_BIKER), color=key)) +
  geom_pointrange(aes(ymin = FCO2_BIKER_low, ymax = FCO2_BIKER_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='grey', method='lm', se=F)+
  geom_line(aes(y=10^(lwr)), color='grey', linetype='dashed', size=1.75) +
  geom_line(aes(y=10^(upr)), color='grey', linetype='dashed', size=1.75) +
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
  annotate("text", label = paste0('RMSE: ', rmse, ' g/m2*dy'), x = 10^-0.4, y = 10^1, size = 5, colour = "black")+
  annotate("text", label = paste0('r2: ', lmr2), x = 10^-0.5, y = 10^0.8, size = 5, colour = "black")

#add some individual river plots (same rivers has randomly sampled from the k600 validation)
#OhioSection3
t <- filter(output, river=='OhioSection3') %>%
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
OhioSection3_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('OhioSection3')

#SacramentoDownstream
t <- filter(output, river=='SacramentoDownstream') %>%
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
SacramentoDownstream_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('SacramentoDownstream')

#SeineDownstream
t <- filter(output, river=='SeineDownstream') %>%
  gather(key=key, value=value, c(FCO2_ulset, FCO2_BIKER))
SeineDownstream_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('SeineDownstream')

#bring it together
legend <- get_legend(OhioSection3_plot + theme(legend.box.margin = margin(0, 3, 3, 5)))
plot_fin <- plot_grid(OhioSection3_plot + theme(legend.position = 'none'), SacramentoDownstream_plot, SeineDownstream_plot, legend, ncol=2,
                      labels=c('b', 'c', 'd', NA))
plot_fin <- plot_grid(flux_plot, plot_fin, ncol=2, labels=c('a', NA))
ggsave('cache/FCO2/FCO2_plot.jpg', plot_fin, width=14, height=8)

#save stats to file--------------------------------------
results_all_rivs <- data.frame('rmse'=rmse, 'r2'=lmr2)
write.csv(results_all_rivs, 'cache/FCO2/fco2_stats_all.csv')
write.csv(stats_by_reach, 'cache/FCO2/fco2_stats_by_river.csv')
