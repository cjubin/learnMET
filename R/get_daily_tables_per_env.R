#' Obtain daily climate data for an environment based on longitude, latitude, 
#' planting and harvest date characterizing this environment.
#' 
#' @param environment \code{character} Name of the environment for which climate
#' data should be extracted.
#'
#' @param info_environments \code{data.frame} with at least 4 columns.
#' \enumerate{
#'   \item year \code{numeric} Year label of the environment
#'   \item location \code{character} Name of the location
#'   \item longitude \code{numeric} longitude of the environment
#'   \item latitude \code{numeric} latitude of the environment
#'   \item planting.date \code{Date} YYYY-MM-DD
#'   \item harvest.date \code{Date} YYYY-MM-DD 
#'   \item IDenv \code{character} ID environment (combination Year x Location) 
#'   \cr
#'   \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location used in the analyses for four 
#'   years, 4 rows should be present (same information with only the value in 
#'   column year changing).} \
#'}
#' 
#' @return a data.frame \code{data.frame} with the following columns extracted 
#' from POWER data, according to requested parameters:
#' \enumerate{
#'   \item LON \code{numeric} 
#'   \item LAT \code{numeric} 
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
#'   \item WS2M \code{numeric} 
#'   \item IDenv \code{character} ID environment for which weather data were
#'   downloaded.
#'   \item length.gs \code{difftime} length in days of the growing season 
#'   for the environment.
#'   }
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export



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
