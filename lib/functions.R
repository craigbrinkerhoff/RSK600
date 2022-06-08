##############
#some functions I will need
##############

# Sum of squares function-------------------------------------------------------
sumsq <- function(x) sum(x^2)

#Schmidt number for CO2------------------------------------------
Sc_co2_func <- function(t){1742 + (-91.24*t) + (2.208*t^2) + (-0.0219*t^3)} #Raymond etal 2012

#Henry's law constant as a function of water temperature [mol/L*atm]
henry_func <- function(temp) {
  A <- 108.3865 #constant
  B <- 0.01985076 #constant
  C <- -6919.53 #constant
  D <- -40.4515 #constant
  E <- 669365 #constant

  temp_k <- temp + 273.15 #[kelvin]
  output <- 10^(A + B*temp_k + C/temp_k + D*log10(temp_k) + E/temp_k^2)

  return(output)
}

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

#Calculate 'observed' k600 via the hydraulically-wide chainsaw model (see '~src\swot_k_model.R' for its validation)
k600_model <- function(depth, slope, vel) {
  return(exp(3.89 + 0.4320*log(g*colMeans(slope, na.rm=T)) + 0.3084*log(colMeans(vel, na.rm=T)) + 0.5226*log(colMeans(depth, na.rm=T))))
  #return(62.82*(g*colMeans(slope, na.rm=T))^(7/16)*colMeans(vel,na.rm=T)^(1/4)*colMeans(depth,na.rm=T)^(9/16)) #REYNOLDS EXTENSION MODEL FOR ED
}
