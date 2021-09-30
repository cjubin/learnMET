#' Quality control on daily weather data
#'
#' @description
#' This function checks range of values for \code{METData} and implements
#' various test on daily weather data (persistence tests, internal
#' consistency tests) provided by the user.
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
#'   weather variable names.
#'   Column names of weather variables must be given as following:
#'   (\strong{For envs with raw daily weather data, T2M, T2M_MIN, T2M_MAX and 
#'   PRECTOTCORR are mandatory to provide and must be given without missing values, 
#'   which implies that any imputation step should be performed before providing
#'   this dataset to the package. 
#'   .}):
#'    \enumerate{
#'     \item T2M \code{numeric} Daily mean temperature (°C)
#'     \item T2M_MIN \code{numeric} Daily minimum temperature (°C)
#'     \item T2M_MAX \code{numeric} Daily maximum temperature (°C)
#'     \item PRECTOTCORR \code{numeric} Daily total precipitation (mm)
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
#'   \strong{
#'   Warning messages are also thrown if some observations do not pass either
#'   the range test, persistence test or the internal consistency test. A
#'   data.frame with dubious values replaced by a "flagged" character is
#'   provided, along with explanations why these data were flagged (column
#'   "reason") is provided as output. None of the flagged values is assigned as
#'   missing values or transformed; therefore we strongly recommend the user to
#'   have a second look at the daily weather data provided and to correct
#'   potential dubious values indicated by the output of the present function.}
#'   \cr
#'   \strong{
#'   Solar radiation or wind data are automatically retrieved from NASA, if they
#'   are not provided without any missing data by the user. As for any other 
#'   weather variable used in this function, these data cannot be only partially
#'   provided (no missing values accepted).}
#'
#' @references
#' \insertRef{zotarelli2010step}{learnMET}
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
qc_raw_weather_data <-
  function(daily_weather_data,
           info_environments,
           path_flagged_values) {
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
        'DOY',
        "T2M",
        "T2M_MIN",
        "T2M_MAX",
        "PRECTOTCORR"
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
        "PRECTOTCORR",
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
    
    # Check which IDenv are provided with raw weather data
    
    envs_with_daily_wdata <- names(table(daily_weather_data$IDenv))
    
    # Check that the dates provided by the user as raw weather data correspond to
    # those from the growing season of the environment.
    
    for (j in envs_with_daily_wdata) {
      int <-
        lubridate::interval(info_environments[info_environments$IDenv == j, 'planting.date'], info_environments[info_environments$IDenv ==
                                                                                                                  j, 'harvest.date'])
      if (!all(daily_weather_data[daily_weather_data$IDenv == j, "YYYYMMDD"] %within%
               int)) {
        stop(
          paste0(
            "The range of dates provided in the raw weather data for ",
            j,
            " does not correspond to the interval of dates between planting and harvest dates."
          )
        )
      }
      
    }
    
    # Create a data.frame to indicate flagged values and the reason why these
    # were flagged.
    
    flagged_values <- daily_weather_data
    flagged_values$reason <- NA
    
    # Check no missing values are present for the main weather variables:
    
    checkmate::assertDataFrame(daily_weather_data[daily_weather_data$IDenv %in% envs_with_daily_wdata, c("T2M", "T2M_MIN", "T2M_MAX", 'PRECTOTCORR')], any.missing = FALSE)
    
    #### QC on precipitation ####
    
    if ('PRECTOTCORR' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$PRECTOTCORR)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$PRECTOTCORR > 500))) {
        warning("Some daily precipitation data sup. to 500 mm, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$PRECTOTCORR < 0))) {
        warning("Some daily precipitation data inf. to 0 mm, which is abnormal.")
      }
      
      flagged_values$PRECTOTCORR[which(flagged_values$PRECTOTCORR > 500)] <-
        'flagged'
      flagged_values$PRECTOTCORR[which(flagged_values$PRECTOTCORR < 0)] <-
        'flagged'
      flagged_values$reason[which(flagged_values$PRECTOTCORR > 500)] <-
        'range_test'
      flagged_values$reason[which(flagged_values$PRECTOTCORR < 0)] <-
        'range_test'
      
      
      
      
    }
    
    
    #### QC on temperature data ####
    
    if ('T2M_MIN' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$T2M_MIN)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$T2M_MIN < (-50)))) {
        warning("Some min. daily temp. inf. to -50, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$T2M_MIN > 30))) {
        warning("Some min. daily temp. sup. to 30, which is abnormal.")
      }
      
      daily_weather_data$T2M_MIN[which(daily_weather_data$T2M_MIN > 30)] <-
        NA
      daily_weather_data$T2M_MIN[which(daily_weather_data$T2M_MIN <
                                         (-50))] <- NA
      
      flagged_values$T2M_MIN[which(flagged_values$T2M_MIN > 30)] <-
        'flagged'
      flagged_values$T2M_MIN[which(flagged_values$T2M_MIN < (-50))] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M_MIN > 30)] <-
        'range_test'
      flagged_values$reason[which(flagged_values$T2M_MIN < (-50))] <-
        'range_test'
      
    }
    
    if ('T2M_MAX' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$T2M_MAX)
      
      if (any(na.omit(daily_weather_data$T2M_MAX < (-40)))) {
        warning("Some max. daily temp. inf. to -50, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$T2M_MAX > 50))) {
        warning("Some max. daily temp. sup. to 30, which is abnormal.")
      }
      
      
      flagged_values$T2M_MAX[which(flagged_values$T2M_MAX < (-40))] <-
        'flagged'
      flagged_values$T2M_MAX[which(flagged_values$T2M_MAX > 50)] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M_MAX < (-40))] <-
        'range_test'
      flagged_values$reason[which(flagged_values$T2M_MAX > 50)] <-
        'range_test'
      
      
    }
    
    if ('T2M' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$T2M_MIN)
      
      if (any(na.omit(daily_weather_data$T2M < (-50)))) {
        warning("Some mean. daily temp. inf. to -50, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$T2M > 50))) {
        warning("Some mean. daily temp. sup. to 50, which is abnormal.")
      }
      
      
      flagged_values$T2M[which(flagged_values$T2M < (-50))] <-
        'flagged'
      flagged_values$T2M[which(flagged_values$T2M > 50)] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M < (-50))] <-
        'range_test'
      flagged_values$reason[which(flagged_values$T2M > 50)] <-
        'range_test'
      
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
        
        flagged_values$T2M_MAX[which(
          daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value
        )] <- 'consistency_test'
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
        
        flagged_values$T2M_MIN[which(
          daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value
        )] <- 'consistency_test'
        
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
        
        flagged_values$T2M[which(
          daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
            daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
            daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value
        )] <- 'persistence_test'
        
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
        
        flagged_values$T2M_MIN[which(
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value
        )] <- 'persistence_test'
        
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
        
        flagged_values$T2M_MAX[which(
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value
        )] <- 'persistence_test'
        
      }
    }
    
    
    #### QC on daily solar radiation data ####
    
    if ('daily_solar_radiation' %in% names(daily_weather_data) & all(!is.na(daily_weather_data$daily_solar_radiation))) {
      checkmate::assert_numeric(daily_weather_data$daily_solar_radiation)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$daily_solar_radiation > 35))) {
        warning("Daily solar radiation exhibits values sup. to 35 MJ/m2/day, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$daily_solar_radiation < 1))) {
        warning("Daily solar radiation exhibits values inf. to 1 MJ/m2/day, which is abnormal.")
      }
      
      flagged_values$daily_solar_radiation[which(daily_weather_data$daily_solar_radiation > 35)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$daily_solar_radiation > 35)] <-
        'range_test'
      
      flagged_values$daily_solar_radiation[which(daily_weather_data$daily_solar_radiation < 1)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$daily_solar_radiation < 1)] <-
        'range_test'
      
      
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
        
        flagged_values$daily_solar_radiation[which(
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
            daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
            daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value
        )] <- 'persistence_test'
        
      }
      
    } else {
      solar_data <- lapply(
        envs_with_daily_wdata,
        FUN = function(x,
                       ...) {
          get_solar_radiation(environment = x,
                              info_environments = info_environments,
                              ...)
        }
      )
    }
    #### QC on relative humidity ####
    
    # 1) Range test
    if ('RH2M' %in% names(daily_weather_data)) {
      
      
      checkmate::assert_numeric(daily_weather_data$RH2M)
      
      
      if (any(na.omit(daily_weather_data$RH2M > 100))) {
        warning("Some relative humidity data sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$RH2M < 0))) {
        warning("Some relative humidity data inf. to 100, which is abnormal.")
      }
      
      flagged_values$RH2M[which(daily_weather_data$RH2M > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M > 100)] <-
        'range_test'
      
      flagged_values$RH2M[which(daily_weather_data$RH2M < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M < 0)] <-
        'range_test'
    }
    
    if ('RH2M_MIN' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$RH2M_MIN)
      
      
      if (any(na.omit(daily_weather_data$RH2M_MIN > 100))) {
        warning("Some min. relative humidity data sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$RH2M_MIN < 0))) {
        warning("Some min. relative humidity data inf. to 100, which is abnormal.")
      }
      
      flagged_values$RH2M_MIN[which(daily_weather_data$RH2M_MIN > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MIN > 100)] <-
        'range_test'
      
      flagged_values$RH2M_MIN[which(daily_weather_data$RH2M_MIN < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MIN < 0)] <-
        'range_test'
    }
    
    if ('RH2M_MAX' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$RH2M_MAX)
      
      
      if (any(na.omit(daily_weather_data$RH2M_MAX > 100))) {
        warning("Some max. relative humidity data sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$RH2M_MAX < 0))) {
        warning("Some max. relative humidity data inf. to 100, which is abnormal.")
      }
      
      flagged_values$RH2M_MAX[which(daily_weather_data$RH2M_MAX > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MAX > 100)] <-
        'range_test'
      
      flagged_values$RH2M_MAX[which(daily_weather_data$RH2M_MAX < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MAX < 0)] <-
        'range_test'
    }
    
    # 2) Persistence test
    
    
    if ('RH2M' %in% names(daily_weather_data)) {
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
        
        flagged_values$RH2M[which(
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
            daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
            daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value
        )] <- 'persistence_test'
        
        
      }
      
    }
    
    if ('RH2M_MIN' %in% names(daily_weather_data)) {
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmin_previous_day_value = lag(RH2M_MIN, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmin_2_days_before_value = lag(RH2M_MIN, n = 2, order_by = c(IDenv)))
      
      
      if (any(
        na.omit(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value
        )
      )) {
        warning(paste(
          'Daily relative humidity remains exactly constant three days in a row.'
        ))
        
        flagged_values$RH2M_MIN[which(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
            daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
            daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value
        )] <- 'persistence_test'
        
        
      }
      
    }
    
    if ('RH2M_MAX' %in% names(daily_weather_data)) {
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmax_previous_day_value = lag(RH2M_MAX, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmax_2_days_before_value = lag(RH2M_MAX, n = 2, order_by = c(IDenv)))
      
      
      if (any(
        na.omit(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value
        )
      )) {
        warning(paste(
          'Daily relative humidity remains exactly constant three days in a row.'
        ))
        
        flagged_values$RH2M_MAX[which(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
            daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
            daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value
        )] <- 'persistence_test'
        
        
      }
      
    }
    #### QC on wind data ####
    if ('WS2M' %in% names(daily_weather_data) & all(!is.na(daily_weather_data$WS2M))) {
      checkmate::assert_numeric(daily_weather_data$WS2M)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$WS2M > 100))) {
        warning("Some daily wind speed data (m/s) sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$WS2M < 0))) {
        warning("Some daily wind speed data (m/s) inf. to 0, which is abnormal.")
      }
      
      flagged_values$WS2M[which(daily_weather_data$WS2M > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$WS2M > 100)] <-
        'range_test'
      
      flagged_values$WS2M[which(daily_weather_data$WS2M < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$WS2M < 0)] <-
        'range_test'
      
      
    } else {
      wind_data <- lapply(
        envs_with_daily_wdata,
        FUN = function(x,
                       ...) {
          get_wind_data(environment = x,
                              info_environments = info_environments,
                              ...)
        }
      )
    }
    
    ## Calculation of the vapor-pressure deficit: difference between the actual
    ## water vapor pressure and the saturation water pressure at a particular
    ## temperature
    
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
          rhmean = daily_weather_data$RH2M
        )
    } else if (all(c('T2M_MIN', 'T2M_MAX') %in% names(daily_weather_data))) {
      mean_saturation_vapor_pressure <-
        get.es(tmin = daily_weather_data$T2M_MIN, tmax = daily_weather_data$T2M_MAX)
      
    }
    if (exists(actual_vapor_pressure) &
        exists(mean_saturation_vapor_pressure)) {
      daily_w_env$vapr_deficit <-
        mean_saturation_vapor_pressure - actual_vapor_pressure
    }
    
    print("QC on daily weather data is done!")
    
    if ("flagged" %in% flagged_values) {
      saveRDS(
        flagged_values,
        file = paste0(
          path_flagged_values,
          '/flagged_values_raw_weather_data.RDS'
        )
      )
    }
    
    return(daily_weather_data)
  }