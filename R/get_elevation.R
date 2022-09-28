#' Obtain elevation data for each field trial based on longitude and latitude
#' data
#'
#' @description
#' A few steps are required to get the elevation data for a given field trial
#' (there might be a more straightforward manner, though):
#'  \enumerate{
#'     \item Get the country name based on the map.where function from the maps
#'     package
#'     \item Convert the country name to ISO code for countries
#'     \item Use the getData function from raster package using the ISO codes
#'     and the longitude and latitude data to obtain elevation
#'  }
#'
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
#'     \item IDenv: \code{character} ID of the environment (location x year)\cr
#'  }
#'
#'  @return elevation_data\code{data.frame} object with 2 columns. \cr
#'  \enumerate{
#'    \item IDenv
#'    \item elevation
#'  }
#'  @export
#'  @author Cathy C. Westhues

get_elevation <- function(info_environments) {
  
  df <- data.table::data.table(info_environments)
  df <- df[, c('longitude', 'latitude')]
  colnames(df) <- c('lon', 'lat')
  rownames(df) <- info_environments$IDenv
  
  iso_codes <- maps::map.where(x=df$lon,
                               y=df$lat)
  
  country_ids <-
    countrycode::countrycode(iso_codes,
                             origin = 'country.name',
                             destination = 'iso3c')
  
  
  x <- raster::getData('alt', 
                       country = country_ids)
  
  elevation <- raster::extract(x, df, method = "bilinear")
  elevation_data <- as.data.frame(cbind(info_environments$IDenv,
                          elevation))
  colnames(elevation_data) <- c('IDenv', 'elevation')
  
  return(elevation_data)
}