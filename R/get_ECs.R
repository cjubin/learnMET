#' Calculates reference ET0 based on the Penman-Monteith model (FAO-56 Method)
#'
#' @description
#' This function calculates the potential evapotranspiration rate from
#' a reference crop canopy (ET0) in mm/d.
#'
#' For these calculations the
analysis by FAO is followed as laid down in the FAO publication
`Guidelines for computing crop water requirements - FAO Irrigation
and drainage paper 56 <http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents>`_
.


get_ECs <- function(METData,customized_growth_intervals=F,nb_intervals=4,equal_length_stages=T) {


  if (METData$compute_ECs != TRUE) {
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
    stop('Planting and harvest dates needed to impute ECs (format Date, y-m-d)')
  }




  # Obtain daily "AG" community daily weather information for each environment


  res_w_daily_all<- lapply(unique(METData$info_environments$IDenv),FUN = function(x) {
    get_daily_tables_env(
      environment = x,
      info_environments = METData$info_environments
    )
  })


  # According to the choice on the method to derive environmental covariates

  if(customized_growth_intervals==F){
    if (equal_length_stages==T){


      get_average=



    }
  }



}











get_daily_tables_env <- function(environment, info_environments) {


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
