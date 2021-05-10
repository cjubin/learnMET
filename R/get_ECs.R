#' Compute environmental covariates for each environment of the MET dataset
#' based on geographical coordinates and a given time frame.
#'
#' @description This function enables to retrieve daily weather data for each
#' environment and derive environmental covariates over time windows, which
#' can be defined in various ways.
#'
#' @param METData A \code{list} object of class `METData` created by the initial
#'   function of the package create_METData().
#'
#' @param unique_EC_by_geno \code{logical} indicates if the environmental
#'   covariates contained in env_data are also genotype-specific (dependent on
#'   the phenology, for instance) or unique for a complete environment.\cr
#'   Default is `FALSE`.
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
#'
#' @param duration_time_window_days \code{numeric} Number of days spanned by a
#'   time window. Default is 15.
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
#' @param path_daily_weather_tables \code{character} Path of the folder where a
#'   RDS object will be created to save the daily weather tables if saved. (Do
#'   not use a Slash after the name of the last folder.)
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#'
#' @export
#'
#'




get_ECs <-
  function(METData,
           unique_EC_by_geno = F,
           fixed_length_time_windows_across_env = T,
           duration_time_window_days = 15,
           fixed_nb_windows_across_env = F,
           nb_windows_intervals = 5,
           customized_growth_intervals = F,
           save_daily_weather_tables = T,
           path_daily_weather_tables,
           ...) {
    if (!METData$compute_ECs) {
      stop('compute_ECs = T required in METData to compute ECs.\n')
    }
    
    if ((fixed_length_time_windows_across_env  &
         fixed_nb_windows_across_env) |
        (!fixed_length_time_windows_across_env  &
         !fixed_nb_windows_across_env)) {
      stop(
        paste(
          'Either the length of time windows must be fixed across all',
          'environments (with fixed_length_time_windows_across_env = T), or',
          'the total number of windows to use across all environments must be',
          'fixed (with fixed_nb_windows_across_env = T). But both cannot be T',
          'or F at the same time.\n'
        )
      )
    }
    
    if (!METData$compute_ECs) {
      stop(
        paste(
          'No computation of weather-based environmental covariates',
          'required. If ECs should be computed, use compute_ECs=TRUE.\n'
        )
      )
    }
    
    if (is.null(METData$info_environments$longitude) ||
        is.null(METData$info_environments$latitude) ||
        is.na(METData$info_environments$latitude) ||
        is.na(METData$info_environments$longitude)) {
      stop('Longitude and latitude needed to impute ECs.\n')
    }
    
    if (is.null(METData$info_environments$harvest.date) ||
        is.null(METData$info_environments$planting.date) ||
        is.na(METData$info_environments$harvest.date) ||
        is.na(METData$info_environments$planting.date)) {
      stop('Planting and harvest dates needed to impute ECs (format Date, YYYY-MM-DD).\n')
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
    
    cat('Daily weather tables downloaded for each environment!')
    # res_w_daily_all: list containing for each element the daily weather table
    # for the time frame given by the user.
    
    saveRDS(
      res_w_daily_all,
      file.path(
        path_daily_weather_tables,
        "daily_weather_tables_list.RDS"
      )
      
    )
    
    
    # According to the choice on the method to derive environmental covariates:
    
    if (fixed_length_time_windows_across_env) {
      # Each EC is computed over a fixed certain number of days, given by the
      # parameter "duration_time_window_days".
      # The maximum number of time windows (e.g. the total number of ECs)
      # is determined by the shortest growing season across all environments.
      
      number_total_fixed_windows <-
        floor(min(sapply(res_w_daily_all, function(x)
          unique(as.numeric(x[, 'length.gs'])))) / duration_time_window_days)
      
      
      ECs_all_envs <-
        lapply(
          res_w_daily_all,
          FUN = function(x) {
            compute_EC_fixed_length_window(
              table_daily_W = x,
              duration_time_window_days = duration_time_window_days,
              number_total_fixed_windows = number_total_fixed_windows,
              ...
            )
          }
        )
      
      test
      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs$location <-
        stringr::str_split(merged_ECs$IDenv, '_', simplify = T)[, 1]
      merged_ECs$year <-
        stringr::str_split(merged_ECs$IDenv, '_', simplify = T)[, 2]
      merged_ECs <-
        merged_ECs[, c('IDenv', 'year', 'location', colnames(merged_ECs)[colnames(merged_ECs) %notin%
                                                                           c('IDenv', 'year', 'location')])]
      
      cat(
        paste(
          'Environmental covariates derived from the daily weather tables'
          ,
          'with a fixed length of time windows in days across environments!'
        )
      )
      
      
    }
    
    
    if (fixed_nb_windows_across_env) {
      ECs_all_envs <-
        lapply(
          res_w_daily_all,
          FUN = function(x) {
            compute_EC_fixed_number_windows(table_daily_W = x,
                                            nb_windows_intervals = nb_windows_intervals,
                                            ...)
          }
        )
      
      
      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs$location <-
        stringr::str_split(merged_ECs$IDenv, '_', simplify = T)[, 1]
      merged_ECs$year <-
        stringr::str_split(merged_ECs$IDenv, '_', simplify = T)[, 2]
      merged_ECs <-
        merged_ECs[, c('IDenv', 'year', 'location', colnames(merged_ECs)[colnames(merged_ECs) %notin%
                                                                           c('IDenv', 'year', 'location')])]
      
      cat(
        paste(
          'Environmental covariates derived from the daily weather tables'
          ,
          'with a fixed number of time windows across environments!'
        )
      )
    }
    
    
    
    # Add to the table env_data, if this table already contains environmental
    # covariates.
    if (!is.null(METData$env_data) & !unique_EC_by_geno) {
      METData$env_data <-
        merge(merged_ECs, METData$env_data[,-c('year', 'location')], by = 'IDenv')
    }
    
    else{
      METData$env_data <- merged_ECs
    }
    METData$ECs_computed <- TRUE
    
    return(METData)
    
  }
