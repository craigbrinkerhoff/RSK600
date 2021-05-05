##############
#some functions I will need
##############

# Sum of squares function-------------------------------------------------------
sumsq <- function(x) sum(x^2)

#Schmidt number for CO2------------------------------------------
Sc_co2_func <- function(t){1742 + (-91.24*t) + (2.208*t^2) + (-0.0219*t^3)} #Raymond etal 2012

#Henry's law function--------------------------------------------
henrys_law_func <- function(t) {-60.2409+93.4517*(100/(t+273.15))+23.3585*log((t+273.15)/100)} #Weiss 1974

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

#I THINK OUT OF DATE---------------------------------
#k600_craig <- function(w, s, v) {
#  widthRegime <- ifelse(mean(w, na.rm = T) < 10 & mean(s) < 0.05, 1,
#                        ifelse(mean(w, na.rm = T) < 10 & mean(s, na.rm=T) >= 0.05, 5,
#                               ifelse(mean(w, na.rm = T) < 50, 2,
#                                      ifelse(mean(w, na.rm = T)< 100, 3, 4))))
#  if (widthRegime == 1){
#    k <- (111.58121 * (g*s * (v))^0.6488131)
#  }
#  if(widthRegime ==2 ){
#    k <- (109.04977 * (g*s * (v))^0.6624354)
#  }
#  if(widthRegime == 3){
#    k <- (31.84344 * (g*s * (v))^0.4366114)
#  }
#  if(widthRegime == 4) {
#    k <- (14.16939 * (g*s * (v))^0.2724834)
#  }
#  if(widthRegime == 5) {
#    k <- (792.63149 * (g*s * (v))^1.3160065)
#  }
#  return(k)
#}

#k600 function
k600_craig <- function(s, v) {
  return(85.10025*(g*s*v)^0.59957)
}

#ko2 function
ko2_craig <- function(depth, slope) {
  Ustar <- sqrt(g*depth*slope)
  return(48*Ustar)
}
