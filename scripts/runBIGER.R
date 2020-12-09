#Creator: Craig Brinkerhoff
#Date: Fall 2020
#Description: run BIGER algorithm on pepsi rivers

library(ncdf4)
library(tidyverse)
library(BIGER)

#model settings-------------------------------
g <- 9.8 #m/s2
mannings_uncertainity <- 0.25 #see monte_carlo_analysis.R for actual uncertainty
  #For these tests, should assume only error is from manning's equation (0.21), NOT k600 model b/c we are assuming this model is 'true'
  #In actual implementation, we want posterior uncertainity to reflect both (which is ~0.30 per our Monte Carlo Analysis)

#set working directory----------------------------------------------------------
setwd('/home/cbrinkerhoff_umass_edu/RS_evasion/BIGER/')

#BAM-stle functions------------------------------------------------------
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
files = list.files('/home/cbrinkerhoff_umass_edu/RS_evasion/pepsi_rivers/10day/', pattern="*.nc", full.names = TRUE)

#run model-------------------------------------------------------
output <- data.frame('river'= NA, 'time'=NA, 'kobs'=NA, 'kest_mean'=NA, 'kest_low'=NA, 'kest_high'=NA, 'kest_sd'=NA)
for (currPepsi in files[1:length(files)]){
  name <- substr(currPepsi, 55, nchar(currPepsi))
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
  dA_obs <- calcdA_mat(W_obs,H_obs) #m2
  
  #calculate observed k ------------------------------------------------------------------------
  V_obs <- Q_obs / area #m/s
  
  #my rule-based regression model for k600 [m/dy]------------------------------------------
  widthRegime <- ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) < 0.05, 1,
                        ifelse(mean(W_obs, na.rm = T) < 10 & mean(S_obs, na.rm=T) >= 0.05, 5,
                            ifelse(mean(W_obs, na.rm = T) < 50, 2,
                               ifelse(mean(W_obs, na.rm = T)< 100, 3, 4))))
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
    
  k_obs <- colMeans(k_obs)

  #run BIGEER------------------------------------------
  data <- biger_data(w=W_obs, s=S_obs, dA=dA_obs)
  priors <- biger_priors(data)
  priors$river_type_priors$sigma_post = matrix(mannings_uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want manning's uncertainity. Real implementation would use full model uncertainty
  kest <- biger_estimate(bigerdata = data, bigerpriors = priors)
  #kest <- biger_k600pred(run_bigee)
  
  #write to file
  temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma)
  output <- rbind(output, temp)
  print(paste0(name, ' done'))
}

output <- output[-1,]
write.csv(output, 'results_10day.csv')
