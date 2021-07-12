######################
#Description: Runs the BIKER algorithm on SWOT-simulated rivers
#Creator: Craig Brinkerhoff
#Date: Winter 2021
#####################

print('running BIKER...')

#function for parallelizing BIKER runs------------------------------
run_BIKER <- function(currPepsi, errFlag) {
  #run model-------------------------------------------------------
  name <- substr(currPepsi, 40, nchar(currPepsi))
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

  #calculate observed k (depending on measurements error)------------------------------------------------------------------------
  if (errFlag == 1){
    Wtrue=ncvar_get(data_in,'Reach_Timeseries/Wtrue') #observed W and S with no measurement error
    Strue = ncvar_get(data_in, 'Reach_Timeseries/Strue')

    #prep pepsi rivers------------------------------------------------
    Strue[Strue<=0]=NA
    Strue[is.na(Strue)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
    Wtrue[Wtrue<0]=NA

    #Some NA handling in slopes
    if (any(apply(Strue,2,sum,na.rm=TRUE) ==0)){
      remove_index =  which((apply(Strue,2,sum,na.rm=TRUE) ==0) ==TRUE)

      Wtrue=Wtrue[,-remove_index]
      Strue=Strue[,-remove_index]
    }

    if (any(apply(Strue,1,sum,na.rm=TRUE) ==0)){
      remove_index =  which((apply(Strue,1,sum,na.rm=TRUE) ==0) ==TRUE)

      Wtrue=Wtrue[-remove_index,]
      Strue=Strue[-remove_index,]
    }

    D_obs <- area / Wtrue
    k_obs <- k600_ustar_craig(D_obs, Strue) #k600 equation
    k_obs <- colMeans(k_obs, na.rm=T)
  }

  else { #no measurement error
    D_obs <- area / W_obs
    k_obs <- k600_ustar_craig(D_obs, S_obs) #k600 equation
    k_obs <- colMeans(k_obs, na.rm=T)
  }


  #Calculate dA matrix from RS W and H-----------------------------------------------------------
  dA_obs <- calcdA_mat(W_obs,H_obs) #[m2]

  #run BIKER------------------------------------------
  data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs)
  priors <- biker_priors(data)
  priors$sigma_model$sigma_post = matrix(uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want Rh uncertainty. Real implementation would use full model uncertainty
  kest <- biker_estimate(bikerdata = data, bikerpriors = priors, meas_err=F,iter = 2000L) #meas err needs to be removed

  #write to file
  if(errFlag == 1){
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=1, 'Wobs'=colMeans(W_obs, na.rm=T), 'Sobs'=colMeans(S_obs, na.rm=T), 'Dobs'=colMeans(area, na.rm=T)/colMeans(W_obs, na.rm=T))
    temp$Rhobs <- (temp$Wobs * temp$Dobs)/(temp$Wobs + 2*temp$Dobs)
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '_err.csv'))
  }
  else{
    temp <- data.frame('river'= rep(name, length(k_obs)), 'time'=kest$time, 'kobs'=k_obs, 'kest_mean'=kest$mean, 'kest_low'=kest$conf.low, 'kest_high'=kest$conf.high, 'kest_sd'=kest$sigma, 'errFlag'=0, 'Wobs'=colMeans(W_obs, na.rm=T), 'Sobs'=colMeans(S_obs, na.rm=T), 'Dobs'=colMeans(area, na.rm=T)/colMeans(W_obs, na.rm=T))
    temp$Rhobs <- (temp$Wobs * temp$Dobs)/(temp$Wobs + 2*temp$Dobs)
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '.csv'))
  }
}

#Run parallelized function------------------------------
files <- list.files('data/Frasson_etal_2021/IdealDataxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 2
files2 <- list.files('data/Durand_etal_2016/xxxxxxxxxxxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 1
files <- c(files, files2)

#print(files)
#results <- run_BIKER(files[32], 0)

system.time(
  results <- mclapply(files, run_BIKER, 0, mc.cores=cores) #run no measurement errors
)

files <- list.files('data/Frasson_etal_2021/FullUncertainty', pattern="*.nc", full.names = TRUE)
system.time(
  results <- mclapply(files, run_BIKER, 1, mc.cores=cores) #run with SWOT measurement errors
)

#Make final full results file----------------------------
df_fin <- list.files(path='cache/validation/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows

write.csv(df_fin, 'cache/validation/BIKER_validation_results.csv')
