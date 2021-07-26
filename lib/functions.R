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

#k600 function using craig's model
k600_ustar_craig <- function(Rh, slope) {
  Ustar <- sqrt(g*Rh*slope)
  return(56.0294*Ustar)
}
