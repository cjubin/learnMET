#' Compute environmental covariates for each environment of the MET dataset.
#'
#' @description
#' This function enables to retrieve daily weather data from NASA POWER source
#' for each environment and derive environmental covariates over non-overlapping
#' time windows, which can be defined in various ways by the user.
#' The user can also provide own daily weather data, even for only part of the
#' total number of environments. For the remaining environments, weather data
#' will be retrieved using the NASA POWER query.
#'
#' @param info_environments a \code{data.frame} with the following columns
#'
#'
#'
#' @param fixed_length_time_windows_across_env \code{logical} indicates if the
#'   growing season lengths should be divided in non-overlapping time windows of
#'   fixed lengths (in days) across all environments.
#'   This implies that the total number of time windows, which need to be common
#'   across all environments, is determined by the shortest growing season
#'   included in the MET environments. This further implies that the total
#'   growing season may not be covered by the environmental predictors for
#'   the longest growing seasons.\cr
#'   Default is `TRUE`.
#
#'
#' @param fixed_nb_windows_across_env \code{logical} indicates if the
#'   growing season lengths should be divided in a fixed number of
#'   non-overlapping windows, which fully cover the growing season of each
#'   environment. This means that these time windows might not be of same length
#'   across environments (but always of same length in one environment). \cr
#'   Default is `FALSE`.
#'
#' @param nb_windows_intervals \code{numeric} Number of time windows covering
#'   the growing season length (common number of time windows across all
#'   environments). Default is 5.
#'
#' @param customized_growth_intervals \code{logical} indicates if the growth
#'   intervals to use are given by user and should be used to define time
#'   windows. \cr
#'   Default is `FALSE`.
#'
#' @param save_daily_weather_tables \code{logical} indicates whether the
#'   daily weather tables should be saved. Default is `TRUE`.
#'
#' @param path_data \code{character} Path of the folder where a
#'   RDS object will be created to save the daily weather tables if saved. (Do
#'   not use a Slash after the name of the last folder.)
#'
#' @param ... Arguments passed to the [compute_EC()] function.
#' \strong{Not all of the aforementioned weather variables need to be
#'    provided. Weather variables which are not provided and needed to compute
#'    environmental covariables will be retrieved from NASA POWER and merged to
#'    the weather data given by the user. If these environmental covariables
#'    should not be used in predictions, the user can specify in the prediction
#'    function the list of environmental covariables to use with the argument
#'    list_env_predictors.}
#'    \strong{data are imputed if missing or assigned to NA after QC}
#' @param method_ECs_intervals GDD
#'
#' @return A \code{data.frame} object containing the weather-based environmental
#'   covariates.
#'
#' @references
#' \insertRef{sparks2018nasapower}{learnMET}
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#'
#' @export
#'
#'




get_ECs <-
  function(info_environments,
           raw_weather_data = NULL,
           method_ECs_intervals = 'fixed_nb_windows_across_env',
           length_minimum_gs = NULL,
           save_daily_weather_tables = T,
           path_data = NULL,
           crop_model = NULL,
           nb_windows_intervals = 10,
           duration_time_window_days = 10,
           base_temperature = 10,
           ...) {
    # Check the path_folder: create if does not exist
    
    
    
    if (!is.null(path_data)) {
      path_data <- file.path(path_data, 'weather_data')
      if (!dir.exists(path_data)) {
        dir.create(path_data, recursive = T)
      }
    }
    
    
    
    if (is.null(info_environments$longitude) ||
        is.null(info_environments$latitude) ||
        is.na(info_environments$latitude) ||
        is.na(info_environments$longitude)) {
      stop('Longitude and latitude needed to impute ECs.\n')
    }
    
    if (is.null(info_environments$harvest.date) ||
        is.null(info_environments$planting.date) ||
        is.na(info_environments$harvest.date) ||
        is.na(info_environments$planting.date)) {
      stop('Planting and harvest dates needed to impute ECs (format Date, YYYY-MM-DD).\n')
    }
    
    # Checking that data are in the past to retrieve weather data
    
    assertive.datetimes::assert_all_are_in_past(x = info_environments_G2F$planting.date)
    assertive.datetimes::assert_all_are_in_past(x = info_environments_G2F$harvest.date)
    
    # Check if raw weather data for some environments are provided.
    # If yes, check which weather variables are provided.
    
    if (!is.null(raw_weather_data)) {
      if (is.null(path_data)) {
        stop(
          'Please indicate a path using path_to_save argument where the results from the QC can be stored.'
        )
      }
      
      print(
        paste(
          'Raw weather data are provided by the user and will be used',
          'to build environmental covariates. If some weather variables required',
          'for computation of ECS are not within the provided dataset, they will',
          'be retrieved and added given the NASA POWER source data.'
        )
      )
      
      raw_weather_data$IDenv <-
        paste0(raw_weather_data$location, '_', raw_weather_data$year)
      raw_weather_data$DOY = as.integer(lubridate::yday(raw_weather_data$YYYYMMDD))
      raw_weather_data <-
        qc_raw_weather_data(
          daily_weather_data = raw_weather_data,
          info_environments = info_environments,
          path_flagged_values = path_data
        )
      
      variables_raw_data <-
        colnames(raw_weather_data)
      list_envs_to_retrieve_all_data <-
        info_environments$IDenv[which(info_environments$IDenv %notin% raw_weather_data$IDenv)]
      
    } else{
      variables_raw_data <- NULL
      list_envs_to_retrieve_all_data <-
        unique(info_environments$IDenv)
    }
    
    
    
    ############################################################################
    # Obtain daily "AG" community daily weather information for each environment
    # using nasapower R package
    ############################################################################
    if (!is.null(list_envs_to_retrieve_all_data)) {
      has_unsuccessful_requests <- TRUE
      counter <- 1
      list_envs_loop <- list_envs_to_retrieve_all_data
      # This is an empty list to which all requested data will be assigned.
      requested_data <-
        vector(mode = "list",
               length = length(list_envs_loop))
      names(requested_data) <- list_envs_loop
      
      # Issues with the NASAPOWER query: it sometimes fail --> use of tryCath and while procedure to ensure weather data for each envrionment
      # are retrieved.
      while (has_unsuccessful_requests) {
        res_w_daily_all <-
          lapply(list_envs_loop,
                 function(environment, ...) {
                   requested_data <- tryCatch({
                     get_daily_tables_per_env(environment = environment,
                                              info_environments = info_environments,
                                              ...)
                   },
                   error = function(e)
                     return(NULL),
                   warning = function(w)
                     return(NULL))
                   
                   requested_data
                 })
        names(res_w_daily_all) <- list_envs_loop
        unsuccessful_request_bool <- vapply(res_w_daily_all,
                                            FUN = is.null,
                                            FUN.VALUE = logical(1))
        
        failed_requests <-
          list_envs_loop[unsuccessful_request_bool]
        good_requests <-
          list_envs_loop[!unsuccessful_request_bool]
        
        list_envs_loop <- list_envs_loop[failed_requests]
        
        
        requested_data[good_requests] <-
          res_w_daily_all[good_requests]
        
        counter <- counter + 1
        
        if (counter == 10) {
          stop("At least one request failed five times.", call. = FALSE)
        }
        
        has_unsuccessful_requests <- any(unsuccessful_request_bool)
      }
      
      
      cat('Daily weather tables downloaded for each environment!\n')
      climate_data_retrieved <- TRUE
      
    } else {
      climate_data_retrieved <- FALSE
    }
    
    #######################################################################
    ## Pre-processing to merge user data + NASA data
    #######################################################################
    
    
    if (is.null(raw_weather_data)) {
      # All weather date were retrieved from NASAPOWER data source
      weather_data_list <- requested_data
    }
    
    if (!is.null(raw_weather_data)) {
      if (length(list_envs_to_retrieve_all_data) > 1) {
        # Missing weather information for some environments is binded to
        # weather information provided by the user for some environments.
        
        raw_weather_data$IDenv <- as.factor(raw_weather_data$IDenv)
        weather_data_list <-
          split(raw_weather_data, raw_weather_data$IDenv)
        weather_data_list <-
          append(weather_data_list, requested_data)
        
        
      }
      
      if (length(list_envs_to_retrieve_all_data) == 0) {
        # All weather date were provided by the user as daily weather data tables.
        
        raw_weather_data$IDenv <- as.factor(raw_weather_data$IDenv)
        weather_data_list <-
          split(raw_weather_data, raw_weather_data$IDenv)
        
      }
    }
    
    # Save daily weather data used to compute ECs
    
    if (save_daily_weather_tables & !is.null(path_data)) {
      saveRDS(weather_data_list,
              file.path(path_data,
                        "daily_weather_tables_list.RDS"))
    }
    #############################################
    # Derivation of EC based on selected method #
    #############################################
    
    
    if (method_ECs_intervals == 'GDD') {
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_gdd(table_daily_W = x,
                           crop_model = crop_model,
                           ...)
          }
        )
      
      
      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs <-
        merged_ECs[, c('IDenv', 'year', 'location', colnames(merged_ECs)[colnames(merged_ECs) %notin%
                                                                           c('IDenv', 'year', 'location')])]
      
      
      
      
    }
    
    if (method_ECs_intervals == 'fixed_length_time_windows_across_env') {
      # Each EC is computed over a fixed certain number of days, given by the
      # parameter "duration_time_window_days".
      # The maximum number of time windows (e.g. the total number of ECs)
      # is determined by the shortest growing season across all environments.
      
      length_minimum_gs <- min(vapply(weather_data_list, function(x)
        unique(as.numeric(x[, 'length.gs'])), numeric(1)))
      
      
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_fixed_length_window(
              table_daily_W = x,
              length_minimum_gs = length_minimum_gs,
              duration_time_window_days = duration_time_window_days,
              base_temperature = base_temperature,
              ...
            )
          }
        )
      
      
      
      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs <-
        merged_ECs[, c('IDenv', 'year', 'location', colnames(merged_ECs)[colnames(merged_ECs) %notin%
                                                                           c('IDenv', 'year', 'location')])]
      
      
      
    }
    
    
    if (method_ECs_intervals == 'fixed_nb_windows_across_env') {
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_fixed_number_windows(
              table_daily_W = x,
              nb_windows_intervals = nb_windows_intervals,
              base_temperature = base_temperature,
              ...
            )
          }
        )
      
      
      merged_ECs <- do.call("rbind", ECs_all_envs)
      
      merged_ECs <-
        merged_ECs[, c('IDenv', 'year', 'location', colnames(merged_ECs)[colnames(merged_ECs) %notin%
                                                                           c('IDenv', 'year', 'location')])]
      
      
    }
    
    
    merged_ECs <- list('ECs' = merged_ECs,
                       'climate_data_retrieved' = climate_data_retrieved)
    
    
    return(merged_ECs)
    
  }
