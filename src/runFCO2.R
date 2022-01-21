###########################
#Description: Uses BIKER and rating curves to calculate FCO2 and median C fluxes from the SWOT simulated rivers
#Craig Brinkerhoff
#Winter 2022
##########################

print('starting FCO2...')

#legend labels------------------
legend_labels <- c('BIKER',
                   'Model using observed \ndepth')

#This model specifically predicts w and d from HG and then obtains velocity by continuity. The model fit for velocity was notably worse.
Brinkerhoff2019_AHG <- read.csv('cache/AHG.csv')

###################
##FUNCTION TO GENERATE FCO2 FROM ALL DEPTH MODELS AND ALL SWOT RIVERS--------------------------------
####################
fun_FCO2_analysis <- function(river) {
  #read in swot rivers
  name <- substr(river, 40, nchar(river))
  name <- substr(name,1,nchar(name)-3)
  data_in =nc_open(river)

  W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
  S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
  area = ncvar_get(data_in, 'Reach_Timeseries/A')
  Q_obs = ncvar_get(data_in, 'Reach_Timeseries/Q')
  reaches <- ncvar_get(data_in, 'River_Info/rch_bnd')

  #Account for models known to not start on 1/1/XXXX.
  start_date <- ncvar_get(data_in, 'Reach_Timeseries/t') #grab first date for river
  start_date <- ifelse(start_date[1] !=1, start_date[1], 1) #if this isn't 1/1 , then use that first provided value. Otherwise, start in janurary.
  river_start_date <- format(as.Date('0000/01/01') + start_date, format='%m') #convert to date formats via lubridate package

  #prep river hydraulic models
  S_obs[S_obs<=0]=NA
  S_obs[is.na(S_obs)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
  W_obs[W_obs<0]=NA
  area[area<0]=NA
  Q_obs[Q_obs<0]=NA

  #Sample every 14 days because CO2 data is biweekly (conveniently, average SWOT overpass frequency is 11 days, so close)
  W_obs = W_obs[,seq(river_start_date, ncol(W_obs), samplingRate)]
  S_obs = S_obs[,seq(river_start_date, ncol(S_obs), samplingRate)]
  Q_obs = Q_obs[,seq(river_start_date, ncol(Q_obs), samplingRate)]
  area = area[,seq(river_start_date, ncol(area), samplingRate)]
  
  if(is.null(dim(W_obs))==1) {
    return(NA)
  }

  #Some NA handling in slopes (borrowed from Frassion etal 2021 BAM runs)-----------------------------------------------------------
  if (any(apply(S_obs,2,sum,na.rm=TRUE) ==0)){ #removes timesteps with all NA slopes
    remove_index =  which((apply(S_obs,2,sum,na.rm=TRUE) ==0) ==TRUE)

    W_obs=W_obs[,-remove_index]
    S_obs=S_obs[,-remove_index]
    Q_obs=Q_obs[,-remove_index]
    area=area[,-remove_index]
  }

  if (any(apply(S_obs,1,sum,na.rm=TRUE) ==0)){ #removes sptial steps with all NA slopes
    remove_index =  which((apply(S_obs,1,sum,na.rm=TRUE) ==0) ==TRUE)

    W_obs=W_obs[-remove_index,]
    S_obs=S_obs[-remove_index,]
    Q_obs=Q_obs[,-remove_index]
    area=area[,-remove_index]
  }

  s <- rbind(s, data.frame('river'=river, meanS=mean(S_obs, na.rm=T)))

  #calculate observed average flow depth and veocity
  V_obs <- Q_obs / area #[m/s]
  D_obs <- area / W_obs #[m]

  #Calculate depth rating curves following suite of rating curve models [m/s]
  D_raymond2012 <- exp(-0.895)*Q_obs^0.294
  D_raymond2013 <- exp(1/(1.86*-1.06))*Q_obs^(1-0.51-0.12)
  d <- matrix(NA, ncol=ncol(D_obs), nrow=nrow(D_obs))
  for(i in 1:nrow(D_obs)) {
    for(k in 1:ncol(D_obs)){
      d[i,k] <- mean(c(D_raymond2012[i,k], D_raymond2013[i,k]))
    }
  }
  D_raymond2013 <- d
  D_brinkerhoff2019 <- Brinkerhoff2019_AHG[2,]$coef*Q_obs^(Brinkerhoff2019_AHG[2,]$slope) #Brinkerhoff etal 2019, estimated in CONUS_CO2 code

  #Calculate velocity rating curves following suite of rating curve models [m/s]
  V_raymond2012 <- exp(-1.64)*Q_obs^0.285
  V_raymond2013 <- exp(-1.06)*V_obs^(0.12)
  v <- matrix(NA, ncol=ncol(V_obs), nrow=nrow(V_obs))
  for(i in 1:nrow(V_obs)) {
    for(k in 1:ncol(V_obs)){
      v[i,k] <- mean(c(V_raymond2012[i,k], V_raymond2013[i,k]))
    }
  }

  V_raymond2013 <- v
  V_brinkerhoff2019 <- (1/(Brinkerhoff2019_AHG[2,]$coef*Brinkerhoff2019_AHG[1,]$coef))*Q_obs^(1-Brinkerhoff2019_AHG[2,]$slope-Brinkerhoff2019_AHG[1,]$slope) #Brinkerhoff etal 2019, estimated in CONUS_CO2 code

  #Generate k600 estimates [m/dy]---------------------------
  k600_obs <- k600_model(D_obs, S_obs, V_obs)
  k600_BIKER <- filter(BIKER_results, river == name)$kest_mean
  k600_BIKER <- k600_BIKER[seq(1, length(k600_BIKER), samplingRate)]
  k600_BIKER_low <- filter(BIKER_results, river == name)$kest_low
  k600_BIKER_low <- k600_BIKER_low[seq(1, length(k600_BIKER_low), samplingRate)]
  k600_BIKER_high <- filter(BIKER_results, river == name)$kest_high
  k600_BIKER_high <- k600_BIKER_high[seq(1, length(k600_BIKER_high), samplingRate)]
  k600_raymond2012 <- k600_model(D_raymond2012, S_obs, V_raymond2012)
  k600_raymond2013 <- k600_model(D_raymond2013, S_obs, V_raymond2013)
  k600_brinkerhoff2019 <- k600_model(D_brinkerhoff2019, S_obs, V_brinkerhoff2019)

  #Calculate FCO2 and Sc using Beauliu data-------------------------------------------
  #sort by starting month here
  Data_sort <- Data[(as.numeric(river_start_date)*2):nrow(Data),]
  co2 <- Data_sort$CO2_uatm[1:length(k600_BIKER)] #measured co2
  Sc_co2 <- Sc_co2_func(Data_sort$Water_temp_C) #Sc from measured water temperature
  Sc_co2 <- Sc_co2[1:length(k600_BIKER)]
  henrys_law <- henry_func(Data_sort$Water_temp_C) #Henry's solubility coefficient from measured water temperature
  henrys_law <- henrys_law[1:length(k600_BIKER)]

  #Convert to kco2 [m/dy]--------------------------
  kco2_obs <- (600/Sc_co2)^(1/2)*k600_obs
  kco2_BIKER <- (600/Sc_co2)^(1/2)*k600_BIKER
  kco2_BIKER_low <- (600/Sc_co2)^(1/2)*k600_BIKER_low
  kco2_BIKER_high <- (600/Sc_co2)^(1/2)*k600_BIKER_high
  kco2_raymond2012 <- (600/Sc_co2)^(1/2)*k600_raymond2012
  kco2_raymond2013 <- (600/Sc_co2)^(1/2)*k600_raymond2013
  kco2_brinkerhoff2019 <- (600/Sc_co2)^(1/2)*k600_brinkerhoff2019

  #Calculate actual CO2 fluxes--------------------------------------------
  FCO2_obs <- kco2_obs*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365)  #[g-C/m2*yr]
  FCO2_BIKER <- kco2_BIKER*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]
  FCO2_BIKER_low <- kco2_BIKER_low*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]
  FCO2_BIKER_high <- kco2_BIKER_high*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]
  FCO2_raymond2012 <- kco2_raymond2012*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]
  FCO2_raymond2013 <- kco2_raymond2013*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]
  FCO2_brinkerhoff2019 <- kco2_brinkerhoff2019*((co2-pCO2_a)*(1/1e6)*henrys_law)*(1/0.001)*molarMass*(365) #[g-C/m2*yr]

  #Calculate river surface areas per subreach and timestep--------------------------------------------
  #first, get reach lengths from the cumulative reach lengths in river data
  r <- reaches[length(reaches)]

  #then calculate a representative surface area
  sa <- mean(W_obs, na.rm = T) * r #[m2]

  #Gather results into single df--------------------------------------------
  for_plot <- data.frame('river' = name, 'timestep'=1:length(FCO2_BIKER), 'sa_m2'=rep(sa, length(FCO2_BIKER)),
                        'FCO2_obs'=FCO2_obs,
                        'FCO2_BIKER'=FCO2_BIKER, 'FCO2_BIKER_low'=FCO2_BIKER_low, 'FCO2_BIKER_high'=FCO2_BIKER_high,
                        'FCO2_raymond2012'=FCO2_raymond2012,
                        'FCO2_raymond2013'=FCO2_raymond2013,
                        'FCO2_brinkerhoff2019'=FCO2_brinkerhoff2019)

    #save results to file--------------------------------------------
    write.csv(for_plot, paste0('cache/FCO2/by_river/results_', name, '.csv'))
}

#...................................................................................................................................

###########
##MAKE SURE BY RIVER RESULTS ARE FOR THIS SESSION--------------
###########
fold <- 'cache/FCO2/by_river'

# get all files in the directories, recursively
f <- list.files(fold, include.dirs = F, full.names = T, recursive = T)
file.remove(f) # remove the files

#create dummy slope df (used in function above)
s <- data.frame('river'=NA, 'meanS'=NA)

#################
##READ IN DATA--------------------------------
#################
#Beaulieu timeseries wrangling
data <- read.csv('data/Beauliu2012.csv', header = TRUE)
data <- data[5:30,] #only keep 26 biweekly sample so they reflect a single year (fall 2009-winter/spring/summer 2008). This amounts to
Data <- select(data, c(Date, CO2, Temperature))
Data$Date <- as.Date(as.character(Data$Date), '%m/%d/%y')
Data$Date <- format(Data$Date, format="%m-%d")
Data <- Data[order(Data$Date),]
colnames(Data) <- c('Date', 'CO2_umol_L', 'Water_temp_C')
Data$Water_temp_C <- as.numeric(as.character(Data$Water_temp_C))

#Beaulieu timeseries unit conversions
Data$CO2_mol_L <- as.numeric(as.character(Data$CO2_umol_L))*0.000001 #mol/L to mol/m3
Data$CO2_uatm <- (Data$CO2_mol_L*1e6) * 1/(henry_func(Data$Water_temp_C)) #mol/m3 to uatm

#(((((quick detour to make beaulieu timeseries plot for supplement)))))))))
beaulieu <- ggplot(data=Data, aes(y=CO2_uatm, x=Date)) +
  geom_point(size=7, color='darkgreen') +
  geom_hline(yintercept = 390, linetype='dashed', size=3)+
  ylim(0,2750)+
  ylab('Water-Side CO2 [uatm]') +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size=17),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave('cache/FCO2/Beaulieu_timeseries.jpg', beaulieu, width=15, height=6)

#BIKER results for the SWOT rivers (only those with no measurement error)--------------------------------------------
BIKER_results <- read.csv('cache/validation/BIKER_validation_results.csv')
BIKER_results <- filter(BIKER_results, errFlag == 0)

#read in SWOT rivers (Frasson etal 2021 + Durand etal 2016), no measurement error--------------------------------------------
files <- list.files('data/Frasson_etal_2021/IdealDataxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 2
files2 <- list.files('data/Durand_etal_2016/xxxxxxxxxxxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 1
files <- c(files, files2)

#########################
##RUN FCO2 ANALYSIS IN PARALLEL------------------------------
##########################
sink("cache/fco2_text_dump.txt") #send all stan outputs to a dump file so they don't muck up the console output
results <- mclapply(files, fun_FCO2_analysis, mc.cores=1)

#debugging option
#results <- fun_FCO2_analysis(i)


#########################
##MAKE FULL RESULTS FILE----------------------------
##########################
df_fin <- list.files(path='cache/FCO2/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows()
write.csv(df_fin, 'cache/FCO2/CO2_results.csv')

sink()
closeAllConnections()
