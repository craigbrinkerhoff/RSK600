#packages--------------------------------------------------------------------
library(BIKER)
library(ncdf4)
library(dplyr)

g <- 9.8
uncertainity <- 0.30

#load in rivers----------------------------------------------------------------
files <- list.files('C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600/data/Frasson_etal_2021/IdealDataxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 2
files2 <- list.files('C:/Users/craig/Documents/OneDrive - University of Massachusetts/Ongoing Projects/RSK600/data/Durand_etal_2016/xxxxxxxxxxxxxxxx', pattern="*.nc", full.names = TRUE) #pepsi 1
files <- c(files, files2)

#grab a dummy river
currPepsi <- files[2]

#extract river name-----------------------------------------------------------
name <- substr(currPepsi, 40, nchar(currPepsi))
name <- substr(name,1,nchar(name)-3)
set.seed(12)
data_in =nc_open(currPepsi)

#extract reach-scale river hydraulics-----------------------------------------------------------
W_obs=ncvar_get(data_in,'Reach_Timeseries/W')
H_obs=ncvar_get(data_in,'Reach_Timeseries/H')
S_obs = ncvar_get(data_in, 'Reach_Timeseries/S')
area = ncvar_get(data_in, 'Reach_Timeseries/A')
Q_obs=ncvar_get(data_in,'Reach_Timeseries/Q')
priorQ <- ncvar_get(data_in, 'River_Info/QWBM')

W_obs=W_obs[1:3,1:4]
H_obs=H_obs[1:3,1:4]
S_obs = S_obs[1:3,1:4]
area = area[1:3,1:4]
Q_obs= Q_obs[1:3,1:4]

#prep river hydraulics-----------------------------------------------------------
S_obs[S_obs<=0]=NA
S_obs[is.na(S_obs)] = 0.000017 #min obs SWOT slope Biancarma etal 2016
W_obs[W_obs<0]=NA
H_obs[H_obs<0]=NA
Q_obs[Q_obs<0]=NA
area[area<0]=NA

#Some NA handling in slopes (borrowed from Frassion etal 2021 BAM runs)-----------------------------------------------------------
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

#Calculate observed k600 with no measurement error
#functions to calculate dA from W and H-----------------------------------------
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

#' #' Calculate partial cross-section area from width and height vectors (time series)
#' #' @param w vector of widths
#' #' @param h vector of heights(FROM MARK)
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

#' #Calculate 'observed' k600 via the hydraulically-wide chainsaw model (see '~src\swot_k_model.R' for its validation)
k600_model <- function(depth, slope, vel) {
  #return(76.4*(g*slope)^(9/16)*(depth)^(11/16)) #CHAINSAW MODEL FOR ES
  return(62.82*(g*slope)^(7/16)*vel^(1/4)*depth^(9/16)) #CHAINSAW MODEL FOR ED
}

 D_obs <- area/W_obs #[m]
 V_obs <- Q_obs/area #[m/s]
 dA_obs <- calcdA_mat(W_obs,H_obs) #[m2]
 k_obs <- k600_model(D_obs, S_obs, V_obs) #k600 equation
 k_obs <- colMeans(k_obs, na.rm=T)

#run BIKER------------------------------------------
data <- biker_data(w=W_obs, s=S_obs, dA=dA_obs, priorQ=as.matrix(priorQ))
priors <- biker_priors(data)
#priors$river_type_priors$logk_hat <- rep(log(6.5), ncol(S_obs))
priors$river_type_priors$logk_sd <- rep(0.30, ncol(W_obs)) #0.748
priors$sigma_model$sigma_post = matrix(uncertainity, nrow=nrow(W_obs), ncol=ncol(W_obs)) #For this validation, we only want Rh uncertainty. Real implementation would use full model uncertainty (calculate in '~src\swot_k_model.R')
posterior <- biker_estimate(bikerdata = data, bikerpriors = priors, meas_err=F,iter = 3000L) #meas err needs to be removed
out <- biker_extract(posterior)

#reconstruct posterior k600
n <- colMeans(matrix(out$n$mean, ncol=ncol(S_obs), nrow=nrow(S_obs), byrow = F))
A0 <- colMeans(matrix(out$A0$mean, ncol=ncol(S_obs), nrow=nrow(S_obs), byrow = F))
reconstructed_k600 <- 61.82*g^(7/16)*colMeans(S_obs, na.rm=T)^(9/16)*(1/n)^(1/4)*((A0+colMeans(dA_obs, na.rm=T))/colMeans(W_obs, na.rm=T))^(35/48)

reconstructed_k600
k_obs
out$k600

