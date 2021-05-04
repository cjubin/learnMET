




get_ECs <-
  function(METData,
           fixed_length_time_windows = T,
           duration_time_window_days = 5,
           same_number_time_windows = F,
           nb_window_intervals = 5,
           customized_growth_intervals = F) {
    
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
      # Each EC is computed over a certain number of days. If the variation of 
      # growing season length 
      
     
    }
    
    
    
  }
