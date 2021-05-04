




get_ECs <-
  function(METData,
           fixed_length_time_windows = T,
           duration_time_window_days = 15,
           unequal_number_days_by_window = F,
           nb_windows_intervals = 5,
           customized_growth_intervals = F,
           base_temperature = 10,
           method_GDD_calculation = 'method_b') {
    
    if ((fixed_length_time_windows & unequal_number_days_by_window) |(!fixed_length_time_windows & !unequal_number_days_by_window )){
      stop('Either the windows must be o')
    }
    
    if (!METData$compute_ECs) {
      stop(
        'No computation of weather-based environmental covariates required. If ECs should be computed, use compute_ECs=TRUE'
      )
    }
    
    if (is.null(METData$info_environments$longitude) ||
        is.null(METData$info_environments$latitude) ||
        is.na(METData$info_environments$latitude) ||
        is.na(METData$info_environments$longitude)) {
      stop('Longitude and latitude needed to impute ECs.')
    }
    
    if (is.null(METData$info_environments$harvest.date) ||
        is.null(METData$info_environments$planting.date) ||
        is.na(METData$info_environments$harvest.date) ||
        is.na(METData$info_environments$planting.date)) {
      stop('Planting and harvest dates needed to impute ECs (format Date, YYYY-MM-DD)')
    }
    
    
    
    
    # Obtain daily "AG" community daily weather information for each environment
    
    
    res_w_daily_all <-
      lapply(
        unique(METData$info_environments$IDenv),
        FUN = function(x) {
          get_daily_tables_per_env(environment = x,
                                   info_environments = METData$info_environments)
        }
      )
    
    
    # According to the choice on the method to derive environmental covariates:
    
    if (fixed_length_time_windows) {
      # Each EC is computed over a certain number of days, given by the parameter
      # "duration_time_window_days".
      # The maximum number of time windows (e.g. the total number of ECs) 
      # is determined by the shortest growing season across all environments. 
  
      number_total_fixed_windows<-  floor(min(sapply(res_w_daily_all, function(x)
        unique(as.numeric(x[, 'length.gs']))))/duration_time_window_days)
    
      
      ECs_all_envs <-
        lapply(
          res_w_daily_all,
          FUN = function(x) {
            compute_EC(table_daily_W = x,
                       duration_time_window_days = duration_time_window_days,
                       base_temperature = base_temperature,
                       method_GDD_calculation = method_GDD_calculation,
                       number_total_fixed_windows = number_total_fixed_windows)
          }
        ) 
     
      
      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs$location <- stringr::str_split(merged_ECs$IDenv,'_',simplify = T)[,1]
      merged_ECs$year <- stringr::str_split(merged_ECs$IDenv,'_',simplify = T)[,2]
      merged_ECs <- merged_ECs[,c('IDenv','year','location',colnames(merged_ECs)[colnames(merged_ECs)%notin%c('IDenv','year','location')])]
      
      
    }
    
    METData$env_data <- merged_ECs
    METData$ECs_computed <- TRUE
    
  }
