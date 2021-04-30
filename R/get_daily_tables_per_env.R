#' Obtain daily climate data for an environment based on longitude, latitude, 
#' planting and harvest date characterizing this environment.
#' 
#' @param environment \code{character} Name of the environment for which climate
#' data should be extracted.
#'
#' @param info_environments \code{data.frame} with at least 4 columns.
#'   First column \code{numeric} with the year label
#'   Second column \code{character} with the location
#'   Third column \code{numeric} with the longitude
#'   Fourth column \code{numeric} with the latitude
#'   
#'   Additional (optional) columns, required if the user wants to download 
#'   weather data with the package (via argument compute_ECs = T):
#'   Fifth column \code{Date} planting.date as YYYY-MM-DD
#'   Sixth column \ode{Date} harvest.date as YYYY-MM-DD
#'   \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location used in the analyses for four 
#'   years, 4 rows should be present (same information with only the value in 
#'   column year changing).}
#'
#' 
#' @return a data.frame \code{data.frame} with the following columns extracted 
#' from POWER data, according to requested parameters:
#'   Column 1 \code{numeric} LON
#'   Column 2 \code{numeric} LAT
#'   Column 3 \code{numeric} YEAR
#'   Column 4 \code{integer} MM
#'   Column 5 \code{integer} DD
#'   Column 6 \code{integer} DOY
#'   Column 7 \code{Date} YYYYMMDD
#'   Column 8 \code{numeric} RH2M
#'   Column 9 \code{numeric} T2M
#'   Column 10 \code{numeric} T2M_MIN
#'   Column 11 \code{numeric} T2M_MAX
#'   Column 12 \code{numeric} PRECTOT
#'   Column 13 \code{numeric} ALLSKY_TOA_SW_DWN
#'   Column 14 \code{numeric} ALLSKY_SFC_SW_DWN
#'   Column 15 \code{numeric} T2MDEW
#'   Column 16 \code{numeric} WS2M
#'   Column 17 IDenv \code{character} ID environment for which weather data were
#'   downloaded.
#'   Column 18 length.gs \code{difftime} length in days of the growing season 
#'   for the environment.



get_daily_tables_per_env <- function(environment, info_environments) {
  
  # Check that the data contain planting and harvest dates
  if (length(info_environments$planting.date)==0){stop("Planting date should be provided")}
  if (length(info_environments$harvest.date)==0){stop("Harvest date should be provided")} 
  
  
  if (!requireNamespace('nasapower', quietly = TRUE)) {utils::install.packages("nasapower")}
  if (!requireNamespace('plyr', quietly = TRUE)) {utils::install.packages("plyr")}
  
  longitude = METData$info_environments[METData$info_environments$IDenv==environment,'longitude']
  latitude = METData$info_environments[METData$info_environments$IDenv==environment,'latitude']
  planting.date = METData$info_environments[METData$info_environments$IDenv==environment,'planting.date']
  harvest.date = METData$info_environments[METData$info_environments$IDenv==environment,'harvest.date']
  length.growing.season = difftime(harvest.date,planting.date,units = 'days')
  
  if(!inherits(planting.date,'Date')||!inherits(harvest.date,'Date')){stop('planting.date and harvest.date should be given as Dates (y-m-d).')}
  
  daily_w_env <- get_power(
    community = "AG",
    lonlat = c(
      longitude,
      latitude
    ),
    pars = c("RH2M", "T2M",'T2M_MIN' ,'T2M_MAX', "PRECTOT","ALLSKY_TOA_SW_DWN", "ALLSKY_SFC_SW_DWN","T2MDEW", "WS2M"),
    dates = c(planting.date,harvest.date) ,
    temporal_average = "DAILY"
  )
  
  daily_w_env[daily_w_env== -99]<-NA
  
  NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
  replace(daily_w_env, TRUE, lapply(daily_w_env, NA2mean))
  
  daily_w_env$IDenv<-environment
  daily_w_env$length.gs<-length.growing.season
  
  return(as.data.frame(daily_w_env))
  
}
