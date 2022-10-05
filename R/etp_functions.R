#' Calculates reference ET0 based on the Penman-Monteith model (FAO-56 Method)
#'
#' @description
#' This function calculates the potential evapotranspiration rate from
#' a reference crop canopy (ET0) in mm/d.
#'
#' For these calculations the
#' analysis by FAO is followed as laid down in the FAO publication
#' `Guidelines for computing crop water requirements - FAO Irrigation
#' and drainage paper 56 <http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents>`_
#' solar_radiation: daily shortwave radiation
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export





penman_monteith_reference_et0 <-
  function(doy = NULL,
           latitude = NULL,
           elevation = NULL,
           tmin = NULL,
           tmax = NULL,
           tmean = NULL,
           solar_radiation = NULL,
           wind_speed = NULL,
           rhmean = NULL,
           rhmax = NULL,
           rhmin = NULL,
           tdew = NULL,
           use_rh = TRUE) {
    # Stefan Boltzmann constant (MJ/m2/d/K4)
    STBC = 4.903E-9
    
    # psychrometric instrument constant (kPa/Celsius)
    psycon = 0.665
    
    # albedo and surface resistance [sec/m] for the reference crop canopy
    REFCFC = 0.23
    CRES = 70
    
    G = 0
    
    # latent heat of evaporation of water [MJ/kg]
    LHVAP = 2.45
    
    # atmospheric pressure (kPa)
    t_kelvin = tmean + 273.16
    patm = 101.3 * ((t_kelvin - (0.0065 * elevation)) / t_kelvin) ** 5.26
    
    # psychrometric constant (kPa/Celsius)
    gamma = psycon * patm * 1.0E-3
    
    # saturated vapour pressure (kPa) at temperature tmean
    sat_vap_pressure_meanT = sat_vap_pressure(tmean)
    
    #the slope of the relationship between saturation vapor pressure and temperature
    delta = (4098 * sat_vap_pressure_meanT) / ((tmean + 237.3) ** 2)
    
    
    # mean saturation vapor pressure for a day computed as the mean between the
    # saturation vapor pressure at the mean daily maximum and
    # minimum air temperatures for that period (kPa)
    sat_vap_pressure_maxT = sat_vap_pressure(tmax)
    sat_vap_pressure_minT = sat_vap_pressure(tmin)
    sat_vap = (sat_vap_pressure_minT + sat_vap_pressure_maxT) / 2
    
    # Actual vapor pressure (ea):
    # Either derived from relative humidity
    # Or defined as the measured vapour pressure (kPa) using tdew
    if (use_rh == TRUE) {
      if (!is.null(rhmin) && !is.null(rhmax)) {
        ea = get.ea(
          tmin = tmin,
          tmax = tmax,
          rhmin = rhmin,
          rhmax = rhmax
        )
      }
      else if ((is.null(rhmin) ||
                is.null(rhmax)) && !is.null(rhmean)) {
        ea = get.ea.with.rhmean(tmin = tmin,
                                tmax = tmax,
                                rhmean = rhmean)
      }
      else{
        esmn <- sat_vap_pressure(tmin)
        ea = esmn
      }
      
    }
    
    if (use_rh == FALSE) {
      if (!is.null(tdew)) {
        tmp = (17.27 * tdew) / (tdew + 237.3)
        ea = 0.6108 * exp(tmp)
      }
      else{
        esmn <- sat_vap_pressure(tmin)
        ea = esmn
      }
    }
    
    # Stefan-Boltzmann law
    stb_max = STBC * ((tmax + 273.16) ** 4)
    stb_min = STBC * ((tmin + 273.16) ** 4)
    # Partial net outgoing long-wave radiation term (MJ/m2/d)
    rnl_tmp = ((stb_min + stb_max) / 2) * (0.34 - 0.14 * sqrt(ea))
    
    clear_sky_radiation = (0.75 + (2e-05 * elevation)) * extraterrestrial_rad(doy = doy,
                                                                              latitude = latitude)
    
    
    rnl = rnl_tmp * (1.35 * (solar_radiation / clear_sky_radiation) - 0.35)
    
    rn = ((1 - REFCFC) * solar_radiation - rnl) / LHVAP
    EA = ((900. / (tmean + 273)) * wind_speed * (sat_vap - ea))
    
    gamma_mod = gamma * (1 + (CRES / 208 * wind_speed))
    
    et0 = (delta * (rn - G)) / (delta + gamma_mod) + (gamma * EA) / (delta +
                                                                       gamma_mod)
    
    et0[et0<0]<-0
    
    
    return(et0)
    
  }



# Extraterrestrial radiation, or Angot radiation

extraterrestrial_rad <- function(doy, latitude) {
  # constants
  rad = 0.0174533
  angle = -4.
  
  
  # declination and solar constant for this day
  dec = -asin(sin(23.45 * rad) * cos(2. * pi * (doy + 10.) / 365.))
  SC  = 1370. * (1. + 0.033 * cos(2. * pi * doy / 365.))
  
  #calculation of daylength from intermediate variables
  # SINLD, COSLD and AOB
  SINLD = sin(rad * latitude) * sin(dec)
  COSLD = cos(rad * latitude) * cos(dec)
  AOB = SINLD / COSLD
  
  
  
  daylh = 12 * (1 + 2 * asin(AOB) / pi)
  DSINB  = 3600. * (daylh * SINLD + 24. * COSLD * sqrt(1. - AOB ** 2) /
                      pi)
  DSINBE = 3600. * (daylh * (SINLD + 0.4 * (SINLD ** 2 + COSLD ** 2 *
                                              0.5)) +
                      12. * COSLD * (2. + 3. * 0.4 * SINLD) * sqrt(1. -
                                                                     AOB ** 2) / pi)
  
  
  
  
  
  # extraterrestrial radiation
  angot.radiation = SC * DSINB * 10E-7
  return(angot.radiation)
  
}
