#' Obtain daily climate data for an environment from NASA POWER data.
#'
#' @description
#' Function downloading daily weather data via the package `nasapower` based on
#' longitude, latitude, planting and harvest date characterizing this
#' environment.
#'
#' @param environment \code{character} Name of the environment for which climate
#' data should be extracted.
#'
#' @param info_environments \code{data.frame} object with at least the 4 first
#'   columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD \cr
#'   }
#'   \strong{Input should be `info_environments`.}
#'   \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location evaluated across four years, 4
#'   rows should be present.}
#'
#' @return a data.frame \code{data.frame} with the following columns extracted
#' from POWER data, according to requested parameters:
#' \enumerate{
#'   \item longitude \code{numeric}
#'   \item latitude \code{numeric}
#'   \item YEAR \code{numeric}
#'   \item MM \code{integer}
#'   \item DD \code{integer}
#'   \item DOY \code{integer}
#'   \item YYYYMMDD \code{Date}
#'   \item RH2M \code{numeric}
#'   \item T2M \code{numeric}
#'   \item T2M_MIN \code{numeric}
#'   \item T2M_MAX \code{numeric}
#'   \item PRECTOT \code{numeric}
#'   \item ALLSKY_TOA_SW_DWN \code{numeric}
#'   \item ALLSKY_SFC_SW_DWN \code{numeric}
#'   \item T2MDEW \code{numeric}
#'   \item IDenv \code{character} ID environment for which weather data were
#'   downloaded.
#'   \item length.gs \code{difftime} length in days of the growing season
#'   for the environment.
#'   }
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export



get_daily_tables_per_env <-
  function(environment,
           info_environments,
           variables_raw_data = NULL,
           which_variables_to_add = NULL,
           ...) {
    
    # Check that the data contain planting and harvest dates
    if (length(info_environments$planting.date) == 0) {
      stop("Planting date should be provided")
    }
    if (length(info_environments$harvest.date) == 0) {
      stop("Harvest date should be provided")
    }
    
    
    if (!requireNamespace('nasapower', quietly = TRUE)) {
      utils::install.packages("nasapower")
    }
    if (!requireNamespace('plyr', quietly = TRUE)) {
      utils::install.packages("plyr")
    }
    
    longitude = info_environments[info_environments$IDenv == environment, 'longitude']
    latitude = info_environments[info_environments$IDenv == environment, 'latitude']
    planting.date = info_environments[info_environments$IDenv == environment, 'planting.date']
    harvest.date = info_environments[info_environments$IDenv == environment, 'harvest.date']
    length.growing.season = difftime(harvest.date, planting.date, units = 'days')
    
    
    list_climatic_variables <-
      c(
        "RH2M",
        "T2M",
        'T2M_MIN' ,
        'T2M_MAX',
        "PRECTOT",
        "ALLSKY_TOA_SW_DWN",
        "ALLSKY_SFC_SW_DWN",
        "T2MDEW"
      )
    
    
    if (!inherits(planting.date, 'Date') ||
        !inherits(harvest.date, 'Date')) {
      stop('planting.date and harvest.date should be given as Dates (y-m-d).')
    }
    
    daily_w_env <- get_power(
      community = "AG",
      lonlat = c(longitude,
                 latitude),
      pars = list_climatic_variables,
      dates = c(planting.date, harvest.date) ,
      temporal_average = "DAILY"
    )
    
    daily_w_env[daily_w_env == -99] <- NA
    daily_w_env$ALLSKY_TOA_SW_DWN[is.na(daily_w_env$ALLSKY_TOA_SW_DWN)] <-
      0
    daily_w_env$ALLSKY_SFC_SW_DWN[is.na(daily_w_env$ALLSKY_SFC_SW_DWN)] <-
      0
    daily_w_env$PRECTOT[is.na(daily_w_env$PRECTOT)] <- 0
    
    NA2mean <- function(x)
      replace(x, is.na(x), mean(x, na.rm = TRUE))
    replace(daily_w_env, TRUE, lapply(daily_w_env, NA2mean))
    
    if (is.null(variables_raw_data)) {
      daily_w_env$vapr <-
        (6.1078 * exp((17.269 * daily_w_env$T2MDEW) / (237.3 + daily_w_env$T2MDEW))) /
        10
    }
    daily_w_env$IDenv <- environment
    daily_w_env$length.gs <- length.growing.season
    colnames(daily_w_env)[which(colnames(daily_w_env) == 'ALLSKY_TOA_SW_DWN')] <-
      "top_atmosphere_insolation"
    colnames(daily_w_env)[which(colnames(daily_w_env) == 'ALLSKY_SFC_SW_DWN')] <-
      "daily_solar_radiation"
    colnames(daily_w_env)[which(colnames(daily_w_env) == 'LON')] <-
      "longitude"
    colnames(daily_w_env)[which(colnames(daily_w_env) == 'LAT')] <-
      "latitude"
    
    
    daily_w_env$location <-
      stringr::str_split(daily_w_env$IDenv, '_', simplify = T)[, 1]
    daily_w_env$year <-
      stringr::str_split(daily_w_env$IDenv, '_', simplify = T)[, 2]
    
    daily_w_env <- arrange(daily_w_env,DOY)
      
    return(as.data.frame(daily_w_env))
    
  }
