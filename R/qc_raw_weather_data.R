#' Quality control on daily weather data
#'
#' @description
#' This function checks range of values for \code{METData} and implements
#' various test on daily weather data (persistence tests, internal
#' consistency tests) provided by the user.
#'
#' @param info_environments \code{data.frame} object with at least the 4 first
#'   columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD
#'     \item elevation: (optional) \code{numeric}
#'     \item IDenv: \code{character} ID of the environment (location x year)\cr
#'   }
#'   \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location evaluated across four years, 4
#'   rows should be present.}
#'
#' @param daily_weather_data a \code{data.frame} which contains the following
#'   mandatory columns:
#'   \enumerate{
#'     \item longitude \code{numeric}
#'     \item latitude \code{numeric}
#'     \item year \code{numeric}
#'     \item location \code{character}
#'     \item YYYYMMDD \code{Date} Date of the daily observation written as
#'       YYYY-MM-DD
#'     \item IDenv \code{character} Environment ID written Location_Year
#'     \item T2M \code{numeric} Average mean temperature (degree Celsius)
#'     \item T2M_MIN \code{numeric} Min. temperature (degree Celsius)
#'     \item T2M_MAX \code{numeric} Max. temperature (degree Celsius)
#'     \item PRECTOTCORR \code{numeric} Total daily precipitation (mm)
#'    }
#'   Additional weather data provided by user must be a subset of the following
#'   weather variable names (= next columns):
#'   (\strong{Any imputation step should be performed before providing
#'   this daily weather dataset to the package. }):
#'    \enumerate{
#'     \item RH2M \code{numeric} Daily mean relative humidity (%)
#'     \item RH2M_MIN \code{numeric} Daily minimum relative humidity (%)
#'     \item RH2M_MAX \code{numeric} Daily maximum relative humidity (%)
#'     \item daily_solar_radiation \code{numeric} daily solar radiation
#'     (MJ/m^2/day)
#'     \item T2MDEW \code{numeric} Dew Point (°C)
#'    }
#'    Default is `NULL`.
#'
#'
#' @param et0 whether evapotranspiration should be calculated. False by default.
#'
#' @param path_flagged_values where to save the file with flagged values to
#'   check on (they are not removed from the data, only indicated in the output
#'   file)
#'
#' @return daily_weather_data a  \code{data.frame} after quality check with the
#'   same columns as before the QC. \cr
#'   Vapor pressure deficit is calculated if T2M_MIN, T2M_MAX, and either
#'   RH2M_MIN + RH2M_MAX  or only RH2M are provided.   \cr
#'   et0 calculated if indicated (et0 = TRUE) . \cr
#'   \strong{
#'   Warning messages are also thrown if some observations do not pass either
#'   the range test, persistence test or the internal consistency test. A
#'   data.frame, with dubious values signaled by a column flagged and with the
#'   corresponding explanation in the column "reason", is provided as output.
#'   None of the flagged values is assigned as missing values or transformed;
#'   therefore we strongly recommend the user to have a second look at the daily
#'   weather data provided and to correct potential dubious values indicated by
#'   the output of the present function.}
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
           path_flagged_values,
           et0 = F) {
    cat("QC on daily weather data starts...\n")
    
    checkmate::assert_data_frame(daily_weather_data, any.missing = FALSE)
    
    checkmate::assert_names(
      colnames(daily_weather_data),
      must.include = c(
        "IDenv",
        "location",
        "year",
        "longitude",
        "latitude",
        "YYYYMMDD",
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
        "MM",
        "DD",
        "YEAR",
        "DOY",
        "RH2M",
        "RH2M_MIN",
        "RH2M_MAX",
        "T2M",
        "T2M_MIN",
        "T2M_MAX",
        "PRECTOTCORR",
        "daily_solar_radiation",
        "sunshine_duration",
        "T2MDEW",
        "WS2M",
        "length.gs",
        "vapr_deficit"
      )
    )
    
    # Check YYYYMMDD class
    checkmate::assert_date(daily_weather_data$YYYYMMDD)
    
    # Get DOY if not provided
    if (!("DOY" %in% colnames(daily_weather_data))) {
      daily_weather_data$DOY <-
        lubridate::yday(daily_weather_data$YYYYMMDD)
    }
    
    daily_weather_data$multiple_obs_per_day <-
      paste0(daily_weather_data$IDenv,
             daily_weather_data$DOY)
    if (duplicated(daily_weather_data$multiple_obs_per_day)) {
      cat(
        "Multiple observations for the same day in the same environment were",
        "found and will be removed to keep 1 obs. per day.\n"
      )
      
      daily_weather_data <-
        daily_weather_data[!duplicated(daily_weather_data$multiple_obs_per_day),]
      
    }
    
    
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
      if (!all(lubridate::`%within%`(daily_weather_data[daily_weather_data$IDenv == j, "YYYYMMDD"],
                                     int))) {
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
    flagged_values$flagged <- NA
    flagged_values$reason <- NA
    
    
    
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
      
      # Flagged values
      flagged_values$flagged[which(flagged_values$PRECTOTCORR > 500)] <-
        'flagged'
      flagged_values$flagged[which(flagged_values$PRECTOTCORR < 0)] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$PRECTOTCORR > 500)] <-
        'range_test_precipitation'
      flagged_values$reason[which(flagged_values$PRECTOTCORR < 0)] <-
        'range_test_precipitation'
      
      
      
      
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
      
      
      flagged_values$flagged[which(flagged_values$T2M_MIN > 30)] <-
        'flagged'
      flagged_values$flagged[which(flagged_values$T2M_MIN < (-50))] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M_MIN > 30)] <-
        'range_test_min_temp'
      flagged_values$reason[which(flagged_values$T2M_MIN < (-50))] <-
        'range_test_min_temp'
      
    }
    
    if ('T2M_MAX' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$T2M_MAX)
      
      if (any(na.omit(daily_weather_data$T2M_MAX < (-40)))) {
        warning("Some max. daily temp. inf. to -40, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$T2M_MAX > 50))) {
        warning("Some max. daily temp. sup. to 50, which is abnormal.")
      }
      
      
      flagged_values$flagged[which(flagged_values$T2M_MAX < (-40))] <-
        'flagged'
      flagged_values$flagged[which(flagged_values$T2M_MAX > 50)] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M_MAX < (-40))] <-
        'range_test_max_temp'
      flagged_values$reason[which(flagged_values$T2M_MAX > 50)] <-
        'range_test_max_temp'
      
      
    }
    
    if ('T2M' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$T2M_MIN)
      
      if (any(na.omit(daily_weather_data$T2M < (-50)))) {
        warning("Some mean. daily temp. inf. to -50, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$T2M > 50))) {
        warning("Some mean. daily temp. sup. to 50, which is abnormal.")
      }
      
      
      flagged_values$flagged[which(flagged_values$T2M < (-50))] <-
        'flagged'
      flagged_values$flagged[which(flagged_values$T2M > 50)] <-
        'flagged'
      
      flagged_values$reason[which(flagged_values$T2M < (-50))] <-
        'range_test_mean_temp'
      flagged_values$reason[which(flagged_values$T2M > 50)] <-
        'range_test_mean_temp'
      
    }
    
    # 2) Internal consistency test
    if (all(c('T2M_MIN', 'T2M_MAX', 'T2M') %in% names(daily_weather_data))) {
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
        dplyr::mutate(., min_previous_day_value = dplyr::lag(T2M_MIN, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., max_previous_day_value = dplyr::lag(T2M_MAX, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      
      if (any(
        na.omit(
          daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(
          paste(
            'Max temperature of day j should be superior to min',
            'temperature of day j-1. Check your data'
          )
        )
        
        flagged_values$flagged[which(
          daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MAX < daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'consistency_test_max_temp'
      }
      if (any(
        na.omit(
          daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(
          paste(
            'Min temperature of day j should be inferior to max',
            'temperature of day j-1. Check your data'
          )
        )
        
        flagged_values$flagged[which(
          daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MIN > daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'consistency_test_min_temp'
        
      }
      
      # 3) Persistence test
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., min_previous_day_value = dplyr::lag(T2M_MIN, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      min_2_days_before_value = dplyr::lag(T2M_MIN, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., max_previous_day_value = dplyr::lag(T2M_MAX, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      max_2_days_before_value = dplyr::lag(T2M_MAX, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., mean_previous_day_value = dplyr::lag(T2M, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., mean_2_days_before_value = dplyr::lag(T2M, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      
      
      if (any(
        na.omit(
          daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
          daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Mean temperature remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
            daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M == daily_weather_data_check$mean_previous_day_value &
            daily_weather_data_check$T2M == daily_weather_data_check$mean_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_mean_temp'
        
      }
      
      
      if (any(
        na.omit(
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Min temperature remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_previous_day_value &
            daily_weather_data_check$T2M_MIN == daily_weather_data_check$min_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_min_temp'
        
      }
      
      
      if (any(
        na.omit(
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Max temperature remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_previous_day_value &
            daily_weather_data_check$T2M_MAX == daily_weather_data_check$max_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_max_temp'
        
      }
    }
    
    
    #### QC on daily solar radiation data ####
    
    if ('daily_solar_radiation' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$daily_solar_radiation, any.missing = F)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$daily_solar_radiation > 35))) {
        warning("Daily solar radiation exhibits values sup. to 35 MJ/m2/day, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$daily_solar_radiation < 1))) {
        warning("Daily solar radiation exhibits values inf. to 1 MJ/m2/day, which is abnormal.")
      }
      
      flagged_values$flagged[which(daily_weather_data$daily_solar_radiation > 35)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$daily_solar_radiation > 35)] <-
        'range_test_solar'
      
      flagged_values$flagged[which(daily_weather_data$daily_solar_radiation < 1)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$daily_solar_radiation < 1)] <-
        'range_test_solar'
      
      
      # 2) Persistence test
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      dsr_previous_day_value = dplyr::lag(daily_solar_radiation, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      dsr_2_days_before_value = dplyr::lag(
                        daily_solar_radiation,
                        n = 2,
                        order_by = c(IDenv)
                      ))
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      if (any(
        na.omit(
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Daily solar radiation remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
            daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_previous_day_value &
            daily_weather_data_check$daily_solar_radiation == daily_weather_data_check$dsr_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_solar'
        
      }
      
    } else {
      cat('Get solar radiation data (weather variable not provided by user)\n')
      solar_data <- lapply(
        envs_with_daily_wdata,
        FUN = function(x,
                       ...) {
          get_solar_radiation(environment = x,
                              info_environments = info_environments,
                              ...)
        }
      )
      
      
      
      has_unsuccessful_requests <- TRUE
      counter <- 1
      list_envs_loop <- envs_with_daily_wdata
      # This is an empty list to which all requested data will be assigned.
      solar_data <-
        vector(mode = "list",
               length = length(list_envs_loop))
      names(solar_data) <- list_envs_loop
      
      # Issues with the NASAPOWER query: it sometimes fail --> use of tryCath and while procedure to ensure weather data for each envrionment
      # are retrieved.
      while (has_unsuccessful_requests) {
        res_w_daily_solar <-
          lapply(list_envs_loop,
                 function(environment, ...) {
                   solar_data <- tryCatch({
                     get_solar_radiation(environment = environment,
                                         info_environments = info_environments,
                                         ...)
                   },
                   error = function(e)
                     return(NULL),
                   warning = function(w)
                     return(NULL))
                   
                   solar_data
                 })
        names(res_w_daily_solar) <- list_envs_loop
        unsuccessful_request_bool <- vapply(res_w_daily_solar,
                                            FUN = is.null,
                                            FUN.VALUE = logical(1))
        
        failed_requests <-
          list_envs_loop[unsuccessful_request_bool]
        good_requests <-
          list_envs_loop[!unsuccessful_request_bool]
        
        list_envs_loop <- list_envs_loop[failed_requests]
        
        
        solar_data[good_requests] <-
          res_w_daily_solar[good_requests]
        
        counter <- counter + 1
        
        if (counter == 15) {
          stop("At least one request failed fifteen times.", call. = FALSE)
        }
        
        has_unsuccessful_requests <- any(unsuccessful_request_bool)
      }
      
      
      solar_data <- as.data.frame(data.table::rbindlist(solar_data))
      daily_weather_data <-
        merge(daily_weather_data, solar_data[, c('IDenv', 'YYYYMMDD', 'daily_solar_radiation')], by = c('IDenv', 'YYYYMMDD'))
      
      
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
      
      flagged_values$flagged[which(daily_weather_data$RH2M > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M > 100)] <-
        'range_test_humidity'
      
      flagged_values$flagged[which(daily_weather_data$RH2M < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M < 0)] <-
        'range_test_humidity'
    }
    
    if ('RH2M_MIN' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$RH2M_MIN)
      
      
      if (any(na.omit(daily_weather_data$RH2M_MIN > 100))) {
        warning("Some min. relative humidity data sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$RH2M_MIN < 0))) {
        warning("Some min. relative humidity data inf. to 100, which is abnormal.")
      }
      
      flagged_values$flagged[which(daily_weather_data$RH2M_MIN > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MIN > 100)] <-
        'range_test_min_humidity'
      
      flagged_values$flagged[which(daily_weather_data$RH2M_MIN < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MIN < 0)] <-
        'range_test_min_humidity'
    }
    
    if ('RH2M_MAX' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$RH2M_MAX)
      
      
      if (any(na.omit(daily_weather_data$RH2M_MAX > 100))) {
        warning("Some max. relative humidity data sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$RH2M_MAX < 0))) {
        warning("Some max. relative humidity data inf. to 100, which is abnormal.")
      }
      
      flagged_values$flagged[which(daily_weather_data$RH2M_MAX > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MAX > 100)] <-
        'range_test_max_humidity'
      
      flagged_values$flagged[which(daily_weather_data$RH2M_MAX < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$RH2M_MAX < 0)] <-
        'range_test_max_humidity'
    }
    
    # 2) Persistence test
    
    
    if ('RH2M' %in% names(daily_weather_data)) {
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rh_previous_day_value = dplyr::lag(RH2M, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rh_2_days_before_value = dplyr::lag(RH2M, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      
      
      if (any(
        na.omit(
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Daily relative humidity remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
            daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M == daily_weather_data_check$rh_previous_day_value &
            daily_weather_data_check$RH2M == daily_weather_data_check$rh_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_mean_humidity'
        
        
      }
      
    }
    
    if ('RH2M_MIN' %in% names(daily_weather_data)) {
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmin_previous_day_value = dplyr::lag(RH2M_MIN, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      rhmin_2_days_before_value = dplyr::lag(RH2M_MIN, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      
      
      
      if (any(
        na.omit(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Daily relative humidity remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
            daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_previous_day_value &
            daily_weather_data_check$RH2M_MIN == daily_weather_data_check$rhmin_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_min_humidity'
        
        
      }
      
    }
    
    if ('RH2M_MAX' %in% names(daily_weather_data)) {
      daily_weather_data_check <- daily_weather_data
      
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., rhmax_previous_day_value = dplyr::lag(RH2M_MAX, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(.,
                      rhmax_2_days_before_value = dplyr::lag(RH2M_MAX, n = 2, order_by = c(IDenv)))
      daily_weather_data_check <-
        daily_weather_data_check %>%
        dplyr::mutate(., IDenv_previous_day = dplyr::lag(IDenv))
      
      
      
      
      if (any(
        na.omit(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value &
          daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )
      )) {
        warning(paste(
          'Daily relative humidity remains exactly constant three days in a row.'
        ))
        
        flagged_values$flagged[which(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
            daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'flagged'
        
        flagged_values$reason[which(
          daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_previous_day_value &
            daily_weather_data_check$RH2M_MAX == daily_weather_data_check$rhmax_2_days_before_value &
            daily_weather_data_check$IDenv == daily_weather_data_check$IDenv_previous_day
        )] <- 'persistence_test_max_humidity'
        
        
      }
      
    }
    #### QC on wind data ####
    if ('WS2M' %in% names(daily_weather_data)) {
      checkmate::assert_numeric(daily_weather_data$WS2M, any.missing = F)
      
      # 1) Range test
      if (any(na.omit(daily_weather_data$WS2M > 100))) {
        warning("Some daily wind speed data (m/s) sup. to 100, which is abnormal.")
      }
      if (any(na.omit(daily_weather_data$WS2M < 0))) {
        warning("Some daily wind speed data (m/s) inf. to 0, which is abnormal.")
      }
      
      flagged_values$flagged[which(daily_weather_data$WS2M > 100)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$WS2M > 100)] <-
        'range_test_wind'
      
      flagged_values$flagged[which(daily_weather_data$WS2M < 0)] <-
        'flagged'
      
      flagged_values$reason[which(daily_weather_data$WS2M < 0)] <-
        'range_test_wind'
      
      
    }
    
    ## Calculation of the vapor-pressure deficit: difference between the actual
    ## water vapor pressure and the saturation water pressure at a particular
    ## temperature
    if ('vapr_deficit' %notin% names(daily_weather_data)) {
      print('vapour pressure deficit calculated from humidity and temp. data')
      if (all(c('T2M_MIN', 'T2M_MAX', "RH2M_MIN", "RH2M_MAX") %in% names(daily_weather_data))) {
        cat(
          paste(
            'Actual vapor pressure (ea) calculated from relative humidity',
            'using RH2M_MIN and RH2M_MAX.\n'
          )
        )
        actual_vapor_pressure <-
          get.ea(
            tmin = daily_weather_data$T2M_MIN,
            tmax = daily_weather_data$T2M_MAX,
            rhmin = daily_weather_data$RH2M_MIN,
            rhmax = daily_weather_data$RH2M_MAX
          )
      } else if (all(c('T2M_MIN', "RH2M_MAX") %in% names(daily_weather_data))) {
        cat(
          paste(
            'Actual vapor pressure (ea) calculated from relative humidity',
            'using only RH2M_MAX.\n'
          )
        )
        actual_vapor_pressure <-
          get.ea.with.rhmax(tmin = daily_weather_data$T2M_MIN,
                            rhmax = daily_weather_data$RH2M_MAX)
      } else if (all(c('T2M_MIN', 'T2M_MAX', "RH2M") %in% names(daily_weather_data))) {
        cat(
          paste(
            'Actual vapor pressure (ea) calculated from relative humidity',
            'using RH2M (mean RH).\n'
          )
        )
        actual_vapor_pressure <-
          get.ea.with.rhmean(
            tmin = daily_weather_data$T2M_MIN,
            tmax = daily_weather_data$T2M_MAX,
            rhmean = daily_weather_data$RH2M
          )
      } else{
        cat(
          paste(
            'Actual vapor pressure (ea) calculated from relative humidity',
            'using T2M_MIN.\n'
          )
        )
        
        actual_vapor_pressure <-
          get.ea.no.RH(tmin = daily_weather_data$T2M_MIN)
      }
      
      mean_saturation_vapor_pressure <-
        get.es(tmin = daily_weather_data$T2M_MIN, tmax = daily_weather_data$T2M_MAX)
      
      
      
      
      daily_weather_data$vapr_deficit <-
        mean_saturation_vapor_pressure - actual_vapor_pressure
    }
    
    
    
    ## Calculation of evapotranspiration
    
    
    
    if (et0) {
      cat('et0 is calculated')
      if ('elevation' %in% colnames(info_environments)) {
        daily_weather_data <-
          plyr::join(daily_weather_data, info_environments[, c('IDenv', 'elevation')], by =
                       'IDenv')
        
      }
      if ('elevation' %notin% colnames(info_environments)) {
        elevation <-
          get_elevation(info_environments = info_environments, path =
                          path_data)[, c('IDenv', 'alt')]
        
        daily_weather_data <-
          plyr::join(daily_weather_data, elevation[, c('IDenv', 'elevation')], by =
                       'IDenv')
        
      }
      if ('RH2M_MAX' %notin% names(daily_weather_data)) {
        daily_weather_data$RH2M_MAX <- NULL
      }
      if ('RH2M_MIN' %notin% names(daily_weather_data)) {
        daily_weather_data$RH2M_MIN <- NULL
      }
      
      daily_weather_data <-
        plyr::join(daily_weather_data, info_environments[, c('IDenv', 'latitude', 'longitude')], by =
                     'IDenv')
      
      daily_weather_data$et0 <-
        penman_monteith_reference_et0(
          doy = daily_weather_data$DOY,
          latitude = daily_weather_data$latitude,
          elevation = daily_weather_data$elevation,
          tmin = daily_weather_data$T2M_MIN,
          tmax = daily_weather_data$T2M_MAX,
          tmean = daily_weather_data$T2M,
          solar_radiation = daily_weather_data$ALLSKY_SFC_SW_DWN ,
          wind_speed = daily_weather_data$WS2M,
          rhmean = daily_weather_data$RH2M,
          rhmax = daily_weather_data$RH2M_MAX,
          rhmin = daily_weather_data$RH2M_MIN,
          tdew = NULL,
          use_rh = TRUE
        )
    }
    
    
    
    cat("QC on daily weather data is done!\n")
    
    
    if (any((flagged_values %>% dplyr::select(-YYYYMMDD)) == 'flagged')) {
      cat('A file with flagged values has been saved in the subfolder weather_data.\n')
      saveRDS(
        flagged_values,
        file = paste0(
          path_flagged_values,
          '/flagged_values_raw_weather_data.RDS'
        )
      )
    }
    
    return(as.data.frame(daily_weather_data))
  }