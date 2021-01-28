#Creator: Craig Brinkerhoff
#Date: Winter 2021
#Description: run BIKER algorithm on pepsi rivers

library(ncdf4)
library(tidyverse)
library(BIKER)

options(scipen=999)

#model settings-------------------------------
g <- 9.8 #m/s2
mannings_uncertainity <- 0.25 #see monte_carlo_analysis.R for actual uncertainty
#For these tests, should assume only error is from manning's equation (0.21), NOT k600 model b/c we are assuming this model is 'true'
#In actual implementation, we want posterior uncertainity to reflect both (which is ~0.30 per our Monte Carlo Analysis)

#set working directory----------------------------------------------------------
setwd('/home/cbrinkerhoff_umass_edu/RS_evasion/BIKER/')

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
#using rivers from Frasson etal 2019 & Rodriquez etal 2020
files = list.files('/home/cbrinkerhoff_umass_edu/RS_evasion/SWOTrivs_Rodriguez_etal_2020', pattern="*.nc", full.names = TRUE)

#run model-------------------------------------------------------
output <- data.frame('river'= NA, 'time'=NA, 'kobs'=NA, 'kest_mean'=NA, 'kest_low'=NA, 'kest_high'=NA, 'kest_sd'=NA)
for (currPepsi in files[2:length(files)]){
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
  
  #Sample every 11 days (swot-style)------------------------
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
  
  k_obs <- colMeans(k_obs, na.rm=T)
  
  #run BIKER------------------------------------------
  data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs)
  priors <- biker_priors(data)
  priors$river_type_priors$sigma_post = matrix(mannings_uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want manning's uncertainity. Real implementation would use full model uncertainty
  
  #Run this to set k600 prior as f(slope) a la Mark
  colSobs <- colMeans(log(S_obs), na.rm=T)
  priors$river_type_priors$logk600_hat = ifelse(colSobs < -4.634, 3.22 + 0.347*colSobs, 6.85 + 1.13*colSobs) #needs to be implemented within BIKER
  priors$river_type_priors$logk600_sd =  rep(1.023, ncol(data$Wobs)) #CV of 100%
  
  #run algorithm
  kest <- biker_estimate(bikerdata = data, bikerpriors = priors)
  
  #write to file
  temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma)
  output <- rbind(output, temp)
  print(paste0(name, ' done'))
}

output <- output[-1,]
write.csv(output, 'results_SWOT_11day_slopeK.csv')
