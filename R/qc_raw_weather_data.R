#' Quality control on daily weather data
#'
#' @description
#' This function checks range of values for \code{METData} and implements
#' various test on the daily weather data (persistence tests, internal
#' consistency tests).
#'
#' @param daily_weather_data a \code{data.frame} which contains the following
#'   mandatory columns:
#'   \enumerate{
#'     \item longitude \code{numeric}
#'     \item latitude \code{numeric}
#'     \item year \code{numeric}
#'     \item location \code{character}
#'     \item YYYYMMDD \code{Date}
#'     \item IDenv \code{character}
#'     \item DOY \code{integer}
#'    }
#'   Available weather data provided by user must be a subset of the following
#'   weather variable names. Colnames must be given as following:
#'    \enumerate{
#'     \item T2M \code{numeric} Daily mean temperature (°C)
#'     \item T2M_MIN \code{numeric} Daily minimum temperature (°C)
#'     \item T2M_MAX \code{numeric} Daily maximum temperature (°C)
#'     \item PRECTOT \code{numeric} Daily total precipitation (mm)
#'     \item RH2M \code{numeric} Daily mean relative humidity (%)
#'     \item RH2M_MIN \code{numeric} Daily minimum relative humidity (%)
#'     \item RH2M_MAX \code{numeric} Daily maximum relative humidity (%)
#'     \item daily_solar_radiation \code{numeric} daily solar radiation
#'     (MJ/m^2/day)
#'     \item top_atmosphere_insolation \code{numeric} Top-of-atmosphere
#'     Insolation (MJ/m^2/day)
#'     \item T2MDEW \code{numeric} Dew Point (°C)
#'    }
#'
#'
#' @return a processed \code{data.frame} after quality check with the same
#'   columns as before the QC. \cr
#'   Vapor pressure deficit is calculated if T2M_MIN, T2M_MAX, and either 
#'   RH2M_MIN + RH2M_MAX  or only RH2M are provided.   \cr
#'   \strong{Values which are outside the range of possible
#'   values are assigned to NA. 
#'   Warning messages are also thrown if some observations do not pass either 
#'   the persistency test or the internal consistency test. In this case, 
#'   concerned values are not flagged nor assigned to NA but we recommend the 
#'   user to have a second look at the daily weather data provided.}
#'   
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
qc_raw_weather_data <- function(daily_weather_data) {
  print("QC on daily weather data starts...")
  
  
  checkmate::assert_names(
    colnames(daily_weather_data),
    must.include = c(
      'IDenv',
      'location',
      'year',
      'longitude',
      'latitude',
      'YYYYMMDD',
      'DOY'
    ),
    subset.of = c(
      "IDenv",
      "location",
      "year",
      "longitude",
      "latitude",
      "YYYYMMDD",
      "DOY",
      "RH2M",
      "T2M",
      "T2M_MIN",
      "T2M_MAX",
      "PRECTOT",
      "top_atmosphere_insolation",
      "daily_solar_radiation",
      "T2MDEW",
      "WS2M"
    )
  )
  
  # Check YYYYMMDD class
  checkmate::assert_date(daily_weather_data$YYYYMMDD)
  
  # Order data.frame
  
  daily_weather_data <-
    dplyr::arrange(daily_weather_data, IDenv, DOY)
  
  # Check which IDenv  
  
  #### QC on precipitation ####
  
  if ('PRECTOT' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$PRECTOT)
    
    # 1) Range test
    if (any(na.omit(daily_weather_data$PRECTOT > 500))) {
      warning(
        "Some daily precipitation data sup. to 500 mm, which is abnormal.
              Data assigned to NA."
      )
    }
    if (any(na.omit(daily_weather_data$PRECTOT < 0))) {
      warning("Some daily precipitation data inf. to 0 mm, which is abnormal.
              Data assigned to NA.")
    }
    
    daily_weather_data$PRECTOT[which(daily_weather_data$PRECTOT > 500)] <-
      NA
    daily_weather_data$PRECTOT[which(daily_weather_data$PRECTOT < 0)] <-
      NA
    
  }
  
  
  #### QC on temperature data ####
  
  if ('T2M_MIN' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$T2M_MIN)
    
    # 1) Range test
    if (any(na.omit(daily_weather_data$T2M_MIN < (-50)))) {
      warning("Some min. daily temp. inf. to -50, which is abnormal.
              Data assigned to NA.")
    }
    if (any(na.omit(daily_weather_data$T2M_MIN > 30))) {
      warning("Some min. daily temp. sup. to 30, which is abnormal.
              Data assigned to NA.")
    }
    
    daily_weather_data$T2M_MIN[which(daily_weather_data$T2M_MIN > 30)] <-
      NA
    daily_weather_data$T2M_MIN[which(daily_weather_data$T2M_MIN <
                                       (-50))] <- NA
    
  }
  
  if ('T2M_MAX' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$T2M_MAX)
    
    if (any(na.omit(daily_weather_data$T2M_MAX < (-40)))) {
      warning("Some max. daily temp. inf. to -50, which is abnormal.
              Data assigned to NA.")
    }
    if (any(na.omit(daily_weather_data$T2M_MAX > 50))) {
      warning("Some max. daily temp. sup. to 30, which is abnormal.
              Data assigned to NA.")
    }
    
    
    daily_weather_data$T2M_MAX[which(daily_weather_data$T2M_MAX > 50)] <-
      NA
    daily_weather_data$T2M_MAX[which(daily_weather_data$T2M_MAX <
                                       (-40))] <- NA
    
  }
  
  if ('T2M' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$T2M_MIN)
    
    if (any(na.omit(daily_weather_data$T2M < (-50)))) {
      warning("Some mean. daily temp. inf. to -50, which is abnormal.
              Data assigned to NA.")
    }
    if (any(na.omit(daily_weather_data$T2M > 50))) {
      warning("Some mean. daily temp. sup. to 50, which is abnormal.
              Data assigned to NA.")
    }
    
    
    daily_weather_data$T2M[which(daily_weather_data$T2M > 50)] <-
      NA
    daily_weather_data$T2M[which(daily_weather_data$T2M <
                                   (-50))] <- NA
    
  }
  
  # 2) Internal consistency test
  if (all(c('T2M_MIN', 'T2M_MAX', 'T2M') %in% names(daily_weather_data))) {
    # 1) Range test
    if (any(na.omit(daily_weather_data$T2M_MAX < daily_weather_data$T2M_MIN))) {
      warning(paste(
        'Max temperature should be superior to min temperature.',
        'Check data.'
      ))
    }
    if (any(na.omit(daily_weather_data$T2M_MAX < daily_weather_data$T2M))) {
      warning(paste(
        'Max temperature should be superior to mean temperature.',
        'Check data.'
      ))
    }
    if (any(na.omit(daily_weather_data$T2M_MIN > daily_weather_data$T2M))) {
      warning(paste(
        'Min temperature should be inferior to mean temperature.',
        'Check data.'
      ))
    }
    
    
    
    daily_weather_data_check <- daily_weather_data
    
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., min_previous_day_value = lag(T2M_MIN, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., max_previous_day_value = lag(T2M_MAX, order_by = c(IDenv)))
    
    if (any(
      na.omit(
        daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value
      )
    )) {
      warning(
        paste(
          'Max temperature of day j should be superior to min',
          'temperature of day j-1. Check your data'
        )
      )
    }
    if (any(
      na.omit(
        daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value
      )
    )) {
      warning(
        paste(
          'Min temperature of day j should be inferior to max',
          'temperature of day j-1. Check your data'
        )
      )
    }
    
    # 3) Persistence test
    daily_weather_data_check <- daily_weather_data
    
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., min_previous_day_value = lag(T2M_MIN, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., min_2_days_before_value = lag(T2M_MIN, n = 2, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., max_previous_day_value = lag(T2M_MAX, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., max_2_days_before_value = lag(T2M_MAX, n = 2, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., mean_previous_day_value = lag(T2M, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., mean_2_days_before_value = lag(T2M, n = 2, order_by = c(IDenv)))
    
    
    if (any(
      na.omit(
        daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
        daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value
      )
    )) {
      warning(paste(
        'Mean temperature remains exactly constant three days in a row.'
      ))
    }
    
    
    if (any(
      na.omit(
        daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
        daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value
      )
    )) {
      warning(paste(
        'Min temperature remains exactly constant three days in a row.'
      ))
    }
    
    
    if (any(
      na.omit(
        daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
        daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value
      )
    )) {
      warning(paste(
        'Max temperature remains exactly constant three days in a row.'
      ))
    }
  }
  
  
  #### QC on daily solar radiation data ####
  
  if ('daily_solar_radiation' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$daily_solar_radiation)
    
    # 1) Range test
    if (any(na.omit(daily_weather_data$daily_solar_radiation > 35))) {
      warning(
        "Daily solar radiation exhibits values sup. to 35 MJ/m2/day, which is abnormal.
              Data assigned to NA."
      )
    }
    if (any(na.omit(daily_weather_data$daily_solar_radiation < 1))) {
      warning(
        "Daily solar radiation exhibits values inf. to 1 MJ/m2/day, which is abnormal.
              Data assigned to NA."
      )
    }
    daily_weather_data$daily_solar_radiation[which(daily_weather_data$daily_solar_radiation > 35)] <-
      NA
    daily_weather_data$daily_solar_radiation[which(daily_weather_data$daily_solar_radiation < 1)] <-
      NA
    
    
    # 2) Persistence test
    daily_weather_data_check <- daily_weather_data
    
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(.,
                    dsr_previous_day_value = lag(daily_solar_radiation, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(.,
                    dsr_2_days_before_value = lag(
                      daily_solar_radiation,
                      n = 2,
                      order_by = c(IDenv)
                    ))
    
    
    if (any(
      na.omit(
        daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
        daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value
      )
    )) {
      warning(paste(
        'Daily solar radiation remains exactly constant three days in a row.'
      ))
    }
    
  }
  
  #### QC on relative humidity ####
  
  if ('RH2M' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$RH2M)
    
    # 1) Range test
    if (any(na.omit(daily_weather_data$RH2M > 100))) {
      warning("Some relative humidity data sup. to 100, which is abnormal.
              Data assigned to NA.")
    }
    if (any(na.omit(daily_weather_data$RH2M < 0))) {
      warning("Some relative humidity data inf. to 100, which is abnormal.
              Data assigned to NA.")
    }
    
    daily_weather_data$RH2M[which(daily_weather_data$RH2M > 100)] <-
      NA
    daily_weather_data$RH2M[which(daily_weather_data$RH2M < 0)] <-
      NA
    
    
    # 2) Persistence test
    daily_weather_data_check <- daily_weather_data
    
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., rh_previous_day_value = lag(RH2M, order_by = c(IDenv)))
    daily_weather_data_check <-
      daily_weather_data_check %>%
      dplyr::mutate(., rh_2_days_before_value = lag(RH2M, n = 2, order_by = c(IDenv)))
    
    
    if (any(
      na.omit(
        daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
        daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value
      )
    )) {
      warning(paste(
        'Daily relative humidity remains exactly constant three days in a row.'
      ))
    }
    
  }
  
  
  #### QC on wind data ####
  if ('WS2M' %in% names(daily_weather_data)) {
    checkmate::assert_numeric(daily_weather_data$WS2M)
    
    # 1) Range test
    if (any(na.omit(daily_weather_data$WS2M > 100))) {
      warning("Some daily wind speed data (m/s) sup. to 100, which is abnormal.")
    }
    if (any(na.omit(daily_weather_data$WS2M < 0))) {
      warning("Some daily wind speed data (m/s) inf. to 0, which is abnormal.")
    }
    
    daily_weather_data$WS2M[which(daily_weather_data$WS2M > 100)] <-
      NA
    daily_weather_data$WS2M[which(daily_weather_data$WS2M < 0)] <-
      NA
    
  }
  
  ## Calculation of the vapor-pressure deficit: difference between the actual
  ## water vapor pressure and the saturation water pressure at a particular
  ## temperature
  if (all(c('T2M_MIN', 'T2M_MAX') %in% names(daily_weather_data))) {
    mean_saturation_vapor_pressure <-
      get.es(tmin = daily_weather_data$T2M_MIN, tmax = daily_weather_data$T2M_MAX)
  }
  if (all(c('T2M_MIN', 'T2M_MAX', "RH2M_MIN", "RH2M_MAX") %in% names(daily_weather_data))) {
    actual_vapor_pressure <-
      get.ea(
        tmin = daily_weather_data$T2M_MIN,
        tmax = daily_weather_data$T2M_MAX,
        rhmin = daily_weather_data$RH2M_MIN,
        rhmax = daily_weather_data$RH2M_MAX
      )
  } else if (all(c('T2M_MIN', 'T2M_MAX', "RH2M") %in% names(daily_weather_data))) {
    actual_vapor_pressure <-
      get.ea.with.rhmean(
        tmin = daily_weather_data$T2M_MIN,
        tmax = daily_weather_data$T2M_MAX,
        rhmean = daily_weather_data$RH2M)
  }
  if (exists(actual_vapor_pressure) &
      exists(mean_saturation_vapor_pressure)) {
    daily_w_env$vapr_deficit <-
      mean_saturation_vapor_pressure - actual_vapor_pressure
  }
  
  print("QC on daily weather data is done!")
  
  
  return(daily_weather_data)
}