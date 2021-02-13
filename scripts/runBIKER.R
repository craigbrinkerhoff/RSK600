library(ncdf4)
library(tidyverse)
library(BIKER)
library(parallel)
library(readr)

#no scientific notation
options(scipen=999)

#model settings-------------------------------
err <- 1 #set to 1 to corrupt S and dA with Durand et al. layover model (approximately)
g <- 9.8 #m/s2
mannings_uncertainity <- 0.25 #see monte_carlo_analysis.R for actual uncertainty
#For these tests, should assume only error is from manning's equation (0.21), NOT k600 model b/c we are assuming this model is 'true'
#In actual implementation, we want posterior uncertainity to reflect both (which is ~0.30 per our Monte Carlo Analysis)

#dA functions------------------------------------------------------
#' @param w Matrix of widths
#' @param h Matrix of heights(FROM MARK)
calcdA_mat <- function(w, h) {
  stopifnot(all(dim(w) == dim(h)))
  dA <- w
  for (i in 1:nrow(dA)) {
    dA[i, ] <- calcdA_vec(w[i, ], h[i, ])
  }
  
  dA
}

#' Calculate partial cross-section area from width and height vectors (time series)
#' @param w vector of widths
#' @param h vector of heights(FROM MARK)
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

#load in SWOT observables-------------------------------------------
#Frasson etal 2019 & Rodriquez etal 2020
files <- list.files('inputs/Frasson_etal_2019', pattern="*.nc", full.names = TRUE)
files <- files[-1] #remove Arial Khan

#function for parallelzing BIKER runs------------------------------
run_BIKER <- function(currPepsi, errFlag) {
  #run model-------------------------------------------------------
  name <- substr(currPepsi, 26, nchar(currPepsi))
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
  S_obs[is.na(S_obs)] = 0.000001 #min obs SWOT slope
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
  
  #Calculate dA matrix from RS W and H-----------------------------------------------------------
  dA_obs <- calcdA_mat(W_obs,H_obs) #[m2]
  
  #calculate observed k ------------------------------------------------------------------------
  V_obs <- Q_obs / area #[m/s]
  
  #my rule-based regression model for k600 [m/dy]
  widthRegime <- ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) < 0.05, 1,
                        ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) >= 0.05, 5,
                               ifelse(mean(W_obs, na.rm = T) < 50, 2,
                                      ifelse(mean(W_obs, na.rm = T)< 100, 3, 4))))
  #Generate observed k600 here  
  if (widthRegime == 1){
    k_obs <- (111.58121 * (g*S_obs * (V_obs))^0.6488131)
  }
  if(widthRegime ==2 ){
    k_obs <- (109.04977 * (g*S_obs * (V_obs))^0.6624354)
  }
  if(widthRegime == 3){
    k_obs <- (31.84344 * (g*S_obs * (V_obs))^0.4366114)
  }
  if(widthRegime == 4) {
    k_obs <- (14.16939 * (g*S_obs * (V_obs))^0.2724834)
  }
  if(widthRegime == 5) {
    k_obs <- (792.63149 * (g*S_obs * (V_obs))^1.3160065)
  }
  
  k_obs <- colMeans(k_obs, na.rm=T)
  
  if(errFlag == 1){
    #Introduce measurement errors using Durand et al. 2020 model------------------------
    set.seed(100)
    for(i in 1:nrow(S_obs)) {
      for(j in 1:ncol(S_obs)){
        S_obs[i,j] <- rnorm(1,S_obs[i,j], 1.7e-5) #km/km
        S_obs[i,j] <- ifelse(S_obs[i,j] <= 0, 0.000001, S_obs[i,j]) #min obs slope Biancarma 2016
      }
    }
    
    for(i in 1:nrow(dA_obs)) {
      for(j in 1:ncol(dA_obs)){
        dA_obs[i,j] <- rnorm(1,dA_obs[i,j], W_obs[i,j]*sqrt(2)*0.104) #m2
      }
    }  
  }
  
  #run BIKER------------------------------------------
  data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs)  
  priors <- biker_priors(data)
  priors$sigma_model$sigma_post = matrix(mannings_uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want manning's uncertainity. Real implementation would use full model uncertainty
  
  kest <- biker_estimate(bikerdata = data, bikerpriors = priors, meas_error = FALSE)
  
  #write to file
  if(errFlag == 1){
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=1)
    write.csv(temp, paste0('outputs/validation/by_river/results_', name, '_err.csv'))
  }
  else{
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=0)
    write.csv(temp, paste0('outputs/validation/by_river/results_', name, '.csv'))
  }
}

#Run parallelezied function------------------------------
system.time(
  results <- mclapply(files, run_BIKER, 0) #run no measurement errors
)
system.time(
  results <- mclapply(files, run_BIKER, 1) #run with SWOT measurement errors
)

#Make final full results file----------------------------
df_fin <- list.files(path='outputs/validation/by_river', full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows

write.csv(df_fin, 'outputs/validation/BIKER_validation_results.csv')