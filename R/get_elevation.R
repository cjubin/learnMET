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

get_elevation <- function(info_environments,path) {
  df <- data.table::data.table(info_environments)
  df2 <- info_environments[, c('longitude', 'latitude')]
  colnames(df2) <- c('lon', 'lat')
  rownames(df2) <- info_environments$IDenv
  
  iso_codes <- maps::map.where(x=df$longitude,y=df$latitude)
  
  
  country_ids <-
    countrycode::countrycode(iso_codes,
                             origin = 'country.name',
                             destination = 'iso3c')
  
  df_raster <- vector(mode = 'list',length = length(unique(country_ids)))
  for (s in 1:length(unique(country_ids))) {
    df_raster[[s]] <-
      raster::getData('alt', country = unique(country_ids)[s],path = path)
  }
  
  df_raster2 <- unlist(df_raster)
  
  
  elevation_data <- vector(mode = 'list', length = length(df_raster2))
  
  for (j in 1:length(df_raster2)) {
    elevation_data[[j]] <-
      cbind(df2, alt = raster::extract(df_raster2[[j]], df2))
    elevation_data[[j]]$IDenv <- row.names(elevation_data[[j]])
    
  }
  elevation_data <- as.data.frame(bind_rows(elevation_data))
  elevation_data <- elevation_data[complete.cases(elevation_data),]
  
  checkmate::assert_data_frame(elevation_data, any.missing = F, nrows = nrow(info_environments))
  
  return(elevation_data)
}
