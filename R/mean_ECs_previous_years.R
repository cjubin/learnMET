mean_ECs_previous_years <- function(year,
                                location,
                                longitude,
                                latitude,
                                planting.date,
                                harvest.date,
                                nb_years = 5,
                                ...){
  
  table_past_env <- structure(list('year' = (year-nb_years):(year-1),
                                   'location' = location,
                                   'longitude'=longitude,
                                   'latitude'=latitude),class='data.frame')
  
  table_past_env$planting.date <- stringr::str_sub(planting.date) 
  
}