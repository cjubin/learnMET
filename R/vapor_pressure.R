#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Actual vapor pressure (ea) derived from relative humidity calculated
#' with rhmin and rhmax.
#' 
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' 
#' @param rhmin \code{numeric} minimum daily relative humidity, %
#' @param rhmax \code{numeric} maximum daily relative humidity, %
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @param tmax \code{numeric} maximum daily air temperature, °C
#' @return Actual vapor pressure (ea) 
#' @export 
#' 
get.ea <- function(rhmin, rhmax, tmin, tmax){
  esmn <- get.esmn(tmin)
  esmx <- get.esmx(tmax) 
  ea <- ((esmn * rhmax/100) + (esmx * rhmin/100)) / 2
  return(ea)
}

#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Actual vapor pressure (ea) derived from relative humidity calculated
#' with rhmean, if rhmin and rhmax not available or bad data quality.
#' 
#' @param rhmean \code{numeric} mean daily relative humidity, %
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @param tmax \code{numeric} maximum daily air temperature, °C
#' @return Actual vapor pressure (ea) 
#' @export
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
get.ea.with.rhmean <- function(tmin, tmax, rhmean){
  esmn <- get.esmn(tmin)
  esmx <- get.esmx(tmax) 
  ea <- (rhmean/100) * ((esmn + esmx) / 2)
  return(ea)
}

#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Actual vapor pressure (ea) derived from relative humidity calculated
#' with only rhmax.
#' 
#' @param rhmax \code{numeric} maximum daily relative humidity, %
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @return Actual vapor pressure (ea) 
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export 

get.ea.with.rhmax <- function(tmin,rhmax){
  esmn <- get.esmn(tmin)
  ea <- (rhmax/100) * esmn 
  return(ea)
}

#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Actual vapor pressure (ea) derived from relative humidity calculated
#' with  tmin if no data related to humidity are available.
#' 
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @return Actual vapor pressure (ea) 
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
get.ea.no.RH <- function(tmin){
  esmn <- get.esmn(tmin) # other fun same as James suggested
  ea <- esmn
  return(ea)
}

#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Formula to calculate Saturation vapor pressure (es) at the daily minimum air 
#' temperature
#' 
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @return saturation vapor pressure at the daily minimum air temperature
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
get.esmn <- function(tmin){
  esmn <- .6108 * exp((17.27 * tmin) / (tmin + 237.3))
  return(esmn)
}


#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Formula to calculate Saturation vapor pressure (es) at the daily maximum air 
#' temperature
#' @param tmax \code{numeric} maximum daily air temperature, °C 
#' @return saturation vapor pressure at the daily maximum air temperature
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
get.esmx <- function(tmax){
  esmx <- .6108 * exp((17.27 * tmax) / (tmax + 237.3))
  return(esmx)
}

#' Formulas to compute vapour pressure deficit according to available data
#' 
#' @description
#' Formula to calculate mean saturation vapor pressure from esmx and esmn
#' 
#' @param tmax \code{numeric} maximum daily air temperature, °C 
#' @param tmin \code{numeric} minimum daily air temperature, °C 
#' @return mean saturation vapor pressure is calculated as the mean between the 
#'   saturation vapor pressure at both the daily maximum and minimum air 
#'   temperatures.
#' @export
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}

get.es <- function(tmin, tmax){
  esmn <- get.esmn(tmin)
  esmx <- get.esmx(tmax)
  es <- (esmn + esmx)/2
  return(es)
}

#' Formulas to compute saturated vapor pressure deficit 
#' 
#' @description
#' Formula to calculate saturation vapor pressure from temperature
#' 
#' @param temp \code{numeric} air temperature, °C 
#' @return saturated vapor pressure 
#' @export
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}

sat_vap_pressure <- function(temp) {
  sat_vap_pressure = (0.6108 * exp(17.27 * temp / (temp + 237.3)))
}
