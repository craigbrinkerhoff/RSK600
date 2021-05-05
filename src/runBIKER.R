######################
#Description: Runs the BIKER algorithm on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Winter 2021
#####################

print('running BIKER...')

#load in SWOT observables-------------------------------------------
#Frasson etal 2019 & Rodriquez etal 2020
files <- list.files('data/Frasson_etal_2019', pattern="*.nc", full.names = TRUE)
files <- files[-1] #remove Arial Khan

#function for parallelizing BIKER runs------------------------------
run_BIKER <- function(currPepsi, errFlag, kmodel) {
  #run model-------------------------------------------------------
  name <- substr(currPepsi, 24, nchar(currPepsi))
  name <- substr(name,1,nchar(name)-3)
  set.seed(12)
  data_in =nc_open(currPepsi)

  W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
  H_obs=ncvar_get(data_in,'Reach_Timeseries/H')
  S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
  area = ncvar_get(data_in, 'Reach_Timeseries/A')
  Q_obs = ncvar_get(data_in, 'Reach_Timeseries/Q')

  #prep pepsi rivers------------------------------------------------
  S_obs[S_obs<=0]=NA
  S_obs[is.na(S_obs)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
  W_obs[W_obs<0]=NA
  H_obs[H_obs<0]=NA
  area[area<0]=NA
  Q_obs[Q_obs<0]=NA

  #Sample every 11 days (SWOT-style)------------------------
  W_obs = W_obs[,seq(1, ncol(W_obs), 11)]
  S_obs = S_obs[,seq(1, ncol(S_obs), 11)]
  H_obs = H_obs[,seq(1, ncol(H_obs), 11)]
  Q_obs = Q_obs[,seq(1, ncol(Q_obs), 11)]
  area = area[,seq(1, ncol(area), 11)]

  #Some NA handling in slopes
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

  #calculate observed k ------------------------------------------------------------------------
  V_obs <- Q_obs / area #[m/s]

  #my rule-based regression model for k600 [m/dy]
#  widthRegime <- ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) < 0.05, 1,
#                        ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) >= 0.05, 5,
#                               ifelse(mean(W_obs, na.rm = T) < 50, 2,
#                                      ifelse(mean(W_obs, na.rm = T)< 100, 3, 4))))
  #Generate observed k600 here
#  if (widthRegime == 1){
#    k_obs <- (111.58121 * (g*S_obs * (V_obs))^0.6488131)
#  }
#  if(widthRegime ==2 ){
#    k_obs <- (109.04977 * (g*S_obs * (V_obs))^0.6624354)
#  }
#  if(widthRegime == 3){
#    k_obs <- (31.84344 * (g*S_obs * (V_obs))^0.4366114)
#  }
#  if(widthRegime == 4) {
#    k_obs <- (14.16939 * (g*S_obs * (V_obs))^0.2724834)
#  }
#  if(widthRegime == 5) {
#    k_obs <- (792.63149 * (g*S_obs * (V_obs))^1.3160065)
#  }

  k_obs <- k600_craig(S_obs, V_obs) #k600 equation
  #D_obs <- area / W_obs
  #k_obs <- ko2_craig(D_obs, S_obs) #k600 equation

  k_obs <- colMeans(k_obs, na.rm=T)

  if(errFlag == 1){
    #Introduce measurement errors using Durand et al. 2020 model------------------------
    set.seed(100)
    for(i in 1:nrow(S_obs)) {
      for(j in 1:ncol(S_obs)){
        S_obs[i,j] <- rnorm(1,S_obs[i,j], 1.7e-5) #Durand etal 2020 [km/km]
        S_obs[i,j] <- ifelse(S_obs[i,j] <= 0, 0.000017, S_obs[i,j]) #min obs slope Biancarma 2016
      }
    }

    for(i in 1:nrow(H_obs)) {
      for(j in 1:ncol(H_obs)){
        H_obs[i,j] <- rnorm(1,H_obs[i,j], 0.104) #m2
      }
    }
  }

  #Calculate dA matrix from RS W and H-----------------------------------------------------------
  dA_obs <- calcdA_mat(W_obs,H_obs) #[m2]

  #run BIKER------------------------------------------
  data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs)
  priors <- biker_priors(data, Kmodel = kmodel)
  priors$sigma_model$sigma_post = matrix(mannings_uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want manning's uncertainty. Real implementation would use full model uncertainty
  kest <- biker_estimate(bikerdata = data, bikerpriors = priors, meas_err=F) #this function needs to be removed because it doesn't work

  #write to file
  if(errFlag == 1){
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=1)
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '_err.csv'))
  }
  else{
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=0)
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '.csv'))
  }
}

#Run parallelized function------------------------------
system.time(
  results <- mclapply(files, run_BIKER, 0, 'k600', mc.cores=cores) #run no measurement errors + k600
)
system.time(
  results <- mclapply(files, run_BIKER, 1, 'k600', mc.cores=cores) #run with SWOT measurement errors + k600
)

#results <- run_BIKER(files[1], 1, 1)

#Make final full results file----------------------------
df_fin <- list.files(path='cache/validation/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows

write.csv(df_fin, 'cache/validation/BIKER_validation_results.csv')
