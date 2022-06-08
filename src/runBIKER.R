######################
#Description: Runs the BIKER algorithm on SWOT-simulated rivers (each river run in parallel)
#Creator: Craig Brinkerhoff
#Date: Winter 2022
#####################

print('running BIKER...')

###############
#FUNCTION THAT RUNS BIKER FOR X RIVER (USED TO RUN ALL SWOT RIVERS IN PARALLEL)
##############
run_BIKER <- function(currPepsi, errFlag) {
  #extract river name from netcdfs
  name <- substr(currPepsi, 40, nchar(currPepsi))
  name <- substr(name,1,nchar(name)-3)
  set.seed(12)
  data_in =nc_open(currPepsi)

  #extract reach-scale river hydraulics from netcdfs
  W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
  H_obs=ncvar_get(data_in,'Reach_Timeseries/H')
  S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
  area = ncvar_get(data_in, 'Reach_Timeseries/A')
  Q_obs=ncvar_get(data_in,'Reach_Timeseries/Q')
  priorQ <- ncvar_get(data_in, 'River_Info/QWBM')

  #prep river hydraulics
  S_obs[S_obs<=0]=NA
  S_obs[is.na(S_obs)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
  W_obs[W_obs<0]=NA
  H_obs[H_obs<0]=NA
  Q_obs[Q_obs<0]=NA
  area[area<0]=NA

  #Some NA handling in slopes (borrowed from Frasson etal 2021 BAM runs)
  if (any(apply(S_obs,2,sum,na.rm=TRUE) ==0)){ #removes timesteps with all NA slopes
    remove_index =  which((apply(S_obs,2,sum,na.rm=TRUE) ==0) ==TRUE)

    W_obs=W_obs[,-remove_index]
    H_obs=H_obs[,-remove_index]
    S_obs=S_obs[,-remove_index]
    Q_obs=Q_obs[,-remove_index]
    area=area[,-remove_index]
  }

  if (any(apply(S_obs,1,sum,na.rm=TRUE) ==0)){ #removes sptial steps with all NA slopes
    remove_index =  which((apply(S_obs,1,sum,na.rm=TRUE) ==0) ==TRUE)

    W_obs=W_obs[-remove_index,]
    H_obs=H_obs[-remove_index,]
    S_obs=S_obs[-remove_index,]
    Q_obs=Q_obs[,-remove_index]
    area=area[,-remove_index]
  }

  #IF measurement error is being run, we need to recalculate observed hydraulics using the 'true' variables in the river netcdfs
      #Frasson et al 2021 SWOT river error model is used here
  if (errFlag == 1){
    Wtrue=ncvar_get(data_in,'Reach_Timeseries/Wtrue') #observed W and S with no measurement error
    Strue = ncvar_get(data_in, 'Reach_Timeseries/Strue')
    Htrue = ncvar_get(data_in, 'Reach_Timeseries/Htrue')

    Strue[Strue<=0]=NA
    Strue[is.na(Strue)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
    Wtrue[Wtrue<0]=NA
    Htrue[Htrue<0]=NA

    #Some NA handling in slopes (borrowed from Frasson etal 2021 BAM runs)
    if (any(apply(Strue,2,sum,na.rm=TRUE) ==0)){ #removes timesteps with all NA slopes
      remove_index =  which((apply(Strue,2,sum,na.rm=TRUE) ==0) ==TRUE)

      Wtrue=Wtrue[,-remove_index]
      Htrue=Htrue[,-remove_index]
      Strue=Strue[,-remove_index]
    }

    if (any(apply(Strue,1,sum,na.rm=TRUE) ==0)){ #removes sptial steps with all NA slopes
      remove_index =  which((apply(Strue,1,sum,na.rm=TRUE) ==0) ==TRUE)

      Wtrue=Wtrue[-remove_index,]
      Htrue=Htrue[,-remove_index]
      Strue=Strue[-remove_index,]
    }

    #calculate observed k600 with no measurement error
    Dtrue <- area/Wtrue
    V_obs <- Q_obs/area
    dA_obs <- calcdA_mat(Wtrue,Htrue) #[m2]

    k_obs <- k600_model(Dtrue, Strue, V_obs) #k600 equation
  }

  #Otherwise, calculate observed k600 with no measurement error
  else {
    D_obs <- area/W_obs #[m]
    V_obs <- Q_obs/area #[m/s]
    dA_obs <- calcdA_mat(W_obs,H_obs) #[m2]

    k_obs <- k600_model(D_obs, S_obs, V_obs) #k600 equation
  }

  #run BIKER
  data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs, priorQ=as.matrix(priorQ))
  priors <- biker_priors(data)
  priors$river_type_priors$logk_sd <- rep(0.25, ncol(W_obs)) #obtained from Bayesian linear regression model
  priors$sigma_model$sigma_post = matrix(uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want Rh uncertainty. Real implementation would use full model uncertainty (calculate in '~src\swot_k_model.R')
  kest <- biker_estimate(bikerdata = data, bikerpriors = priors, meas_err=F,iter = 3000L) #measurement error is currently in devolpment, don't use for this study

  #write results to file
  if(errFlag == 1){
    temp <- data.frame('river'= rep(name, length(k_obs)),
                       'time'=kest$time,
                       'kobs'=k_obs,
                       'kest_mean'=kest$mean,
                       'kest_low'=kest$conf.low,
                       'kest_high'=kest$conf.high,
                       'kest_sd'=kest$sigma,
                       'kprior' = exp(priors$river_type_priors$logk_hat),
                       'errFlag'=1,
                       'Wobs'=colMeans(Wtrue, na.rm=T),
                       'Sobs'=colMeans(Strue, na.rm=T),
                       'Dobs'=colMeans(area, na.rm=T)/colMeans(Wtrue, na.rm=T),
                       'dAobs'=colMeans(dA_obs, na.rm=T))
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '_err.csv'))
  }
  else{
    temp <- data.frame('river'= rep(name, length(k_obs)),
                       'time'=kest$time,
                       'kobs'=k_obs,
                       'kest_mean'=kest$mean,
                       'kest_low'=kest$conf.low,
                       'kest_high'=kest$conf.high,
                       'kest_sd'=kest$sigma,
                       'kprior' = exp(priors$river_type_priors$logk_hat),
                       'errFlag'=0,
                       'Wobs'=colMeans(W_obs, na.rm=T),
                       'Sobs'=colMeans(S_obs, na.rm=T),
                       'Dobs'=colMeans(area, na.rm=T)/colMeans(W_obs, na.rm=T),
                       'dAobs'=colMeans(dA_obs, na.rm=T))
    write.csv(temp, paste0('cache/validation/by_river/results_', name, '.csv'))
  }
}

#...................................................................................................................................

#############
#RUN BIKER IN PARALLEL------------------------------
############
set.seed(143)


#run with no measurement error-----------------
  #(x's are placeholders so that I can grab river names from file paths easily- these means total num of characters is equivalent to the measurement error filepath below)
files <- list.files('data/Frasson_etal_2021/IdealDataxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 2
files2 <- list.files('data/Durand_etal_2016/xxxxxxxxxxxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 1
files <- c(files, files2)

sink("cache/biker_text_dump.txt") #send all stan outputs to a dump file so they don't muck up the console output
results <- mclapply(files, run_BIKER, 0, mc.cores=cores)

files <- list.files('data/Frasson_etal_2021/FullUncertainty', pattern="*.nc", full.names = TRUE) #run with SWOT measurement errors
results <- mclapply(files, run_BIKER, 1, mc.cores=cores)

#debugging option
#results <- run_BIKER(files[32], 0)
#break

#################
#CONCATENATE ALL BY-RIVER RESULTS INTO A SINGLE FILE----------------------------
################
df_fin <- list.files(path='cache/validation/by_river', full.names = TRUE) %>%
  lapply(read_csv) %>%
  bind_rows()
write.csv(df_fin, 'cache/validation/BIKER_validation_results.csv')
sink()
