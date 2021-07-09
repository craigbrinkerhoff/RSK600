###########################
#Description: Uses BIKER and rating curves to calculate FCO2 and bulk C fluxes from the SWOT simulated rivers
#Craig Brinkerhoff
#Winter 2021
##########################

print('starting FCO2...')

#legend labels------------------
legend_labels <- c('BIKER',
                   'Model Using Observed \nshear velocity')

###################
##FUNCTION TO GENERATE FCO2 FROM ALL DEPTH MODELS AND ALL SWOT RIVERS
####################
fun_FCO2_analysis <- function(river) {
  #read in swot rivers
  name <- substr(river, 40, nchar(river))
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
  D_raymond2012 <- exp(-0.895)*Q_obs^0.294

  #Calculate velocity using global rating curve from 1/2 the raymond 2013 model (used in Lauerwald 2015 for river orders 4+)
  V_raymond2013 <- exp(-1.06)*Q_obs^0.12
  D_raymond2013 <- exp(1/(1.86*-1.06))*Q_obs^(1-0.51-0.12)

  V_lauerwald2015 <- V_raymond2013
  D_lauerwald2015 <- D_raymond2013

  #avgerage of the models was the actual implementation in Raymond 2013
  d <- data.frame(V_raymond2012, V_raymond2013)
  V_raymond2013 <- rowMeans(d)

  d <- data.frame(D_raymond2012, D_raymond2013)
  D_raymond2013 <- rowMeans(d)

  #Horgby 2019 model for mountain streams. This should be a bad estimate...
  V_horgby2019 <- 0.668*Q_obs^0.365
  D_horgby2019 <- 0.298*Q_obs^0.222

  #eD
  eD <- g * V_obs * S_obs

  #ustar
  Ustar <- sqrt(g*S_obs*D_obs)

  #Generate k600 estimates---------------------------
  k600_obs <- k600_ustar_craig(D_obs, S_obs)
  k600_obs <- colMeans(k600_obs, na.rm=T)
  k600_BIKER <- filter(BIKER_results, river == name)$kest_mean
  k600_BIKER <- k600_BIKER[seq(1, length(k600_BIKER), 11)]
  k600_BIKER_low <- filter(BIKER_results, river == name)$kest_low
  k600_BIKER_low <- k600_BIKER_low[seq(1, length(k600_BIKER_low), 11)]
  k600_BIKER_high <- filter(BIKER_results, river == name)$kest_high
  k600_BIKER_high <- k600_BIKER_high[seq(1, length(k600_BIKER_high), 11)]
  k600_raymond2012 <- k600_ustar_craig(D_raymond2012, S_obs)
  k600_raymond2012 <- colMeans(k600_raymond2012, na.rm=T)
  k600_raymond2013 <- k600_ustar_craig(D_raymond2013, S_obs)
  k600_raymond2013 <- colMeans(k600_raymond2013, na.rm=T)
  k600_lauerwald2015 <- k600_ustar_craig(D_lauerwald2015, S_obs)
  k600_lauerwald2015 <- colMeans(k600_lauerwald2015, na.rm=T)
  k600_horgby2019 <- k600_ustar_craig(D_horgby2019, S_obs)
  k600_horgby2019 <- colMeans(k600_horgby2019, na.rm=T)

  #Calculate FCO2 and Sc using Beauliu data-------------------------------------------
  #sort by starting mmonth here
  Data_sort <- Data[(as.numeric(river_start_date)*2):nrow(Data),]
  co2 <- Data_sort$CO2_uatm[1:length(k600_BIKER)] #measured co2
  Sc_co2 <- Sc_co2_func(Data_sort$Water_temp_C) #Sc from measured water temperature
  Sc_co2 <- Sc_co2[1:length(k600_BIKER)]
  henrys_law <- henrys_law_func(Data_sort$Water_temp_C) #Henry's solubility coefficient from measured water temperature
  henrys_law <- henrys_law[1:length(k600_BIKER)]

  #Convert to kco2--------------------------
  kco2_obs <- (Sc_co2/600)^(-1/2)/k600_obs
  kco2_BIKER <- (Sc_co2/600)^(-1/2)/k600_BIKER
  kco2_raymond2012 <- (Sc_co2/600)^(-1/2)/k600_raymond2012
  kco2_raymond2013 <- (Sc_co2/600)^(-1/2)/k600_raymond2013
  kco2_lauerwald2015 <- (Sc_co2/600)^(-1/2)/k600_lauerwald2015
  kco2_horgby2019 <- (Sc_co2/600)^(-1/2)/k600_horgby2019
  kco2_BIKER_low <- (Sc_co2/600)^(-1/2)/k600_BIKER_low
  kco2_BIKER_high <- (Sc_co2/600)^(-1/2)/k600_BIKER_high

  #convert k for o2 at 20deg to kco2 at x deg
#  kco2_obs <- (Sc_co2/530)^(-1/2)*k600_obs
#  kco2_BIKER <- (Sc_co2/530)^(-1/2)*k600_BIKER
#  kco2_raymond2012 <- (Sc_co2/530)^(-1/2)*k600_raymond2012
#  kco2_raymond2013 <- (Sc_co2/530)^(-1/2)*k600_raymond2013
#  kco2_lauerwald2015 <- (Sc_co2/530)^(-1/2)*k600_lauerwald2015
#  kco2_horgby2019 <- (Sc_co2/530)^(-1/2)*k600_horgby2019
#  kco2_BIKER_low <- (Sc_co2/530)^(-1/2)*k600_BIKER_low
#  kco2_BIKER_high <- (Sc_co2/530)^(-1/2)*k600_BIKER_high

  #Obtain actual CO2 fluxes--------------------------------------------
  FCO2_obs <- kco2_obs*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_BIKER <- kco2_BIKER*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2012 <- kco2_raymond2012*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_raymond2013 <- kco2_raymond2013*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_lauerwald2015 <- kco2_lauerwald2015*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_horgby2019 <- kco2_horgby2019*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

  FCO2_BIKER_low <- kco2_BIKER_low*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2
  FCO2_BIKER_high <- kco2_BIKER_high*((co2-pCO2_a)*exp(henrys_law)*molarMass*1/1000000*1/0.001) #g/m2*dy includes conversion from uatm to mg/L of CO2

  #Surface area
  #first, get reach lengths from the cumulative reach lengths in river data
  r <- reaches
  for (i in 2:length(reaches)) {
    r[i] <- reaches[i] - reaches[i-1]
  }
  r <- r[-1]

  #then calculate surface area
  sa <- sum(rowMeans(W_obs, na.rm = T) * r) #m2

  #plot------------------------
  for_plot <- data.frame('river' = name, 'timestep'=1:length(FCO2_BIKER), 'sa'=sa,
                        'FCO2_obs'=FCO2_obs,
                        'FCO2_BIKER'=FCO2_BIKER, 'FCO2_BIKER_low'=FCO2_BIKER_low, 'FCO2_BIKER_high'=FCO2_BIKER_high,
                        'FCO2_raymond2012'=FCO2_raymond2012,
                        'FCO2_raymond2013'=FCO2_raymond2013,
                        'FCO2_lauerwald2015'=FCO2_lauerwald2015,
                        'FCO2_horgby2019'=FCO2_horgby2019)

    #write results to file----------------
    write.csv(for_plot, paste0('cache/FCO2/by_river/results_', name, '.csv'))
}

s <- data.frame('river'=NA, 'meanS'=NA)

#################
##READ IN DATA--------------------------------
#################
#Beaulieu 2012 Co2 data-------------------------------------
data <- read.csv('data/Beauliu2012.csv', header = TRUE)
data <- data[5:30,] #only keep 26 biweekly sample so they reflect a single year (fall 2009-winter/spring/summer 2008). This amounts to
Data <- select(data, c(Date, CO2, Temperature))
Data$Date <- as.Date(as.character(Data$Date), '%m/%d/%y')
Data$Date <- format(Data$Date, format="%m-%d")
Data <- Data[order(Data$Date),]
colnames(Data) <- c('Date', 'CO2_umol_L', 'Water_temp_C')
Data$Water_temp_C <- as.numeric(as.character(Data$Water_temp_C))
Data$CO2_umol_L <- as.numeric(as.character(Data$CO2_umol_L))
Data$CO2_uatm <- Data$CO2_umol_L / exp(henrys_law_func(Data$Water_temp_C))

#BIKER results (only those with no measurement error)-----------------------------------------------
BIKER_results <- read.csv('cache/validation/BIKER_validation_results.csv')
BIKER_results <- filter(BIKER_results, errFlag == 0)

#SWOT rivers (Frasson etal 2021 + Durand etal 2016)----------------------------------------------------
files <- list.files('data/Frasson_etal_2021/IdealDataxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 2
files <- files[-1] #remove Arial Khan
files2 <- list.files('data/Durand_etal_2016/xxxxxxxxxxxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 1
files <- c(files, files2)

#########################
##RUN FCO2 ANALYSIS IN PARALLEL------------------------------
##########################
system.time(
  results <- mclapply(files, fun_FCO2_analysis, mc.cores=cores)
)

#results <- fun_FCO2_analysis(files[14])

#########################
##MAKE FULL RESULTS FILE----------------------------
##########################
df_fin <- list.files(path='cache/FCO2/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows
write.csv(df_fin, 'cache/FCO2/CO2_results.csv')

##################
##PLOT BEAULIEU 2012 TIMESERIES--------------------------------
###################
beaulieu <- ggplot(data=Data, aes(y=CO2_uatm, x=Date)) +
  geom_point(size=5, color='darkgreen') +
  ylab('Water-Side CO2 [uatm]')
ggsave('cache/FCO2/Beaulieu_timeseries.jpg', beaulieu, width=15, height=6)

#################
##GENERATE RESULTS AND FIGURES----------------
################
#read in full results
output <- read.csv('cache/FCO2/CO2_results.csv')
output$river <- as.character(output$river)

total_sa <- output[!duplicated(output$sa), ] #m2
total_sa <- sum(total_sa$sa)

#CO2 data is only for 29 samples, so ignore k values beyond that (applies to a few rivers)
output <- drop_na(output)

#Calculate total mass fluxes of CO2 from rivers------------------------------------------------
massFluxes <- gather(output, key=key, value=value, c(FCO2_obs, FCO2_BIKER, FCO2_raymond2012, FCO2_raymond2013, FCO2_lauerwald2015, FCO2_horgby2019)) %>%
  group_by(key) %>%
  summarise(sumFCO2 = (sum(value, na.rm=T)/nrow(output)) * total_sa * 1e-9 * (12.011/44.01) * 365) #gG-C/yr from all rivers
write.csv(massFluxes, 'cache/FCO2/bulkFluxes.csv')

##########################
##PLOT BULK EFFLUX AND BY-RIVER METRICS PLOT----------------------------------------------------------------------
###########################
#bar plots of mass fluxes----------------------------------------
massFluxes$key <- factor(massFluxes$key,levels = c("FCO2_BIKER", "FCO2_horgby2019", "FCO2_lauerwald2015", "FCO2_raymond2012", "FCO2_raymond2013", "FCO2_obs"))
barPlots <- ggplot(massFluxes, aes(y=sumFCO2, x=key, fill=key)) +
  geom_bar(stat='identity', color='black', size=1.2) +
  geom_hline(yintercept = massFluxes[massFluxes$key=='FCO2_obs', ]$sumFCO2, size=1.2, linetype='dashed') +
  ylab('Bulk C Efflux [gG-C/yr]') +
  xlab('Velocity Model') +
  scale_fill_brewer(palette = 'Set1', name='', labels=c('BIKER', 'Horgby \n2019', 'Lauerwald \n2015', 'Raymond \n2012', 'Raymond \n2013', 'Observed')) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#Calculate by-river error metrics-----------------------------------
for_plot <- gather(output, key=key, value=value, c(FCO2_raymond2012, FCO2_raymond2013, FCO2_BIKER, FCO2_lauerwald2015, FCO2_horgby2019))
stats_by_reach <- group_by(for_plot, river, key) %>%
  summarise(kge = hydroGOF::KGE(value, FCO2_obs),
            nrmse = sqrt(mean((FCO2_obs - value)^2)) / mean(FCO2_obs, na.rm=T),
            rBIAS =   mean(value - FCO2_obs) / mean(FCO2_obs, na.rm=T),
            rrmse =   sqrt(mean((value - FCO2_obs)^2 / FCO2_obs^2, na.rm=T)))

#plot by-river error metrics---------------------------------------
plot_stats <- gather(stats_by_reach, key=key2, value=value, c('nrmse', 'rBIAS', 'rrmse', 'kge'))
plotRivs <- ggplot(plot_stats, aes(x=key2, y=value, fill=key)) +
  geom_boxplot(size=1) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_hline(yintercept=-0.41, linetype='dashed') +
  xlab('Metric') +
  ylab('Value') +
  coord_cartesian(ylim = c(-2,1.5))+
  scale_fill_brewer(palette = 'Set1', name='') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        legend.title = element_text(size=17, face='bold'))

#Bring these alllllll together and save-----------------------------
FCO2_models <- plot_grid(barPlots, plotRivs, ncol=2, labels='auto')
ggsave('cache/FCO2/FCO2_models.jpg', FCO2_models, width=13, height=6)

#######################
##PLOT FCO2 VALIDATION ACROSS ALL RIVERS AND TIMESTEPS-------------------------
#######################
#Calculate error stats--------------------------------
lm <- lm(log10(FCO2_BIKER)~log10(FCO2_obs), data=output)
lmr2 <- round(summary(lm((FCO2_BIKER)~(FCO2_obs), data=output))$r.squared,2)
rmse <- round((Metrics::rmse((output$FCO2_obs), (output$FCO2_BIKER))), 2) #g/m2/dy

predInts <- predict(lm, interval='prediction')
output <- cbind(output, predInts)

#plot validation-----------------------------------------
flux_plot<- ggplot(output, aes(x=(FCO2_obs), y=(FCO2_BIKER), color=key)) +
  geom_pointrange(aes(ymin = FCO2_BIKER_low, ymax = FCO2_BIKER_high), fatten=10, fill='#1b9e77', pch=21, color='black') +
  geom_smooth(size=2, color='grey', method='lm', se=F)+
  geom_line(aes(y=10^(lwr)), color='grey', linetype='dashed', size=1.75) +
  geom_line(aes(y=10^(upr)), color='grey', linetype='dashed', size=1.75) +
  geom_abline(size=2, linetype='dashed', color='black') +
  xlab('FCO2 via observed \nshear velocity [g/m2*dy]') +
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
  annotate("text", label = paste0('RMSE: ', rmse, ' g/m2*dy'), x = 10^0, y = 10^1.5, size = 5, colour = "black")+
  annotate("text", label = paste0('r2: ', lmr2), x = 10^0, y = 10^1.2, size = 5, colour = "black")

#Plot some individual river timesereis (same as used in the ko2 validation)-----------------------
#MissouriDownstream----------------
t <- filter(output, river=='MissouriDownstream') %>%
  gather(key=key, value=value, c(FCO2_obs, FCO2_BIKER))
MissouriDownstream_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('MissouriDownstream')

#OhioSection8-------------------
t <- filter(output, river=='OhioSection8') %>%
  gather(key=key, value=value, c(FCO2_obs, FCO2_BIKER))
OhioSection8_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('OhioSection8')

#Connecticut---------------------
t <- filter(output, river=='Connecticut') %>%
  gather(key=key, value=value, c(FCO2_obs, FCO2_BIKER))
Connecticut_plot <- ggplot(t, aes(x=timestep, y=value, color=key)) +
  geom_line(size=1.5) +
  geom_pointrange(aes(ymin=FCO2_BIKER_low, ymax=FCO2_BIKER_high), fatten=10)+
  ylab('FCO2 [g/m2*dy]') +
  xlab('Timestep') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette='Set2', name='', labels=legend_labels) +
  ggtitle('Connecticut')

#bring it alllll together and save-----------------------------------
legend <- get_legend(MissouriDownstream_plot + theme(legend.box.margin = margin(0, 3, 3, 5)))
plot_fin <- plot_grid(MissouriDownstream_plot + theme(legend.position = 'none'), OhioSection8_plot, Connecticut_plot, legend, ncol=2,
                      labels=c('b', 'c', 'd', NA))
plot_fin <- plot_grid(flux_plot, plot_fin, ncol=2, labels=c('a', NA))
ggsave('cache/FCO2/FCO2_plot.jpg', plot_fin, width=14, height=8)

##################
##SAVE STATS TO FILES FOR MS--------------------------------------
#################
results_all_rivs <- data.frame('rmse'=rmse, 'r2'=lmr2)
write.csv(results_all_rivs, 'cache/FCO2/fco2_stats_all.csv')
write.csv(stats_by_reach, 'cache/FCO2/fco2_stats_by_river.csv')
