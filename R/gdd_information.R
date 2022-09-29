#' Internal function of [compute_EC_gdd())]
#' @description
#' Get the table of estimated growth stage based on GDD (in Celsius degrees)
#' Only maize (Zea mays) and hard wheat (Triticum durum) currently implemented.
#' @param crop_model name of the crop model to be used to estimate growth stage
#'   based on growing degree days accumulation. \cr
#'   Current possible values are: 'maizehybrid1700' and 'hardwheatUS',
#' @return a \code{data.frame} with:
#'  \enumerate{
#'     \item Stage: \code{character} NAme of the stage
#'     \item GDD: \code{numeric} Accumulated GDDs to reach this growth stage
#'  }
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export

gdd_information <- function(crop_model){
  if (crop_model == 'maizehybrid1700') {
    data("GDD_maize1700")
    base_temperature = 10
    max_temperature = 30
    return(list(GDD_maize1700,base_temperature,max_temperature))
  }
  if (crop_model == 'maizehybrid1300') {
    data("GDD_maize1300")
    base_temperature = 10
    max_temperature = 30
    return(list(GDD_maize1300,base_temperature,max_temperature))
  }
  if (crop_model == 'maizehybrid1500') {
    data("GDD_maize1500")
    base_temperature = 10
    max_temperature = 30
    return(list(GDD_maize1500,base_temperature,max_temperature))
  }
  if (crop_model == 'maizehybrid1900') {
    data("GDD_maize1900")
    base_temperature = 10
    max_temperature = 30
    return(list(GDD_maize1900,base_temperature,max_temperature))
  }
  if (crop_model == 'maizehybrid1800') {
    data("GDD_maize1800")
    base_temperature = 10
    max_temperature = 30
    return(list(GDD_maize1800,base_temperature,max_temperature))
  }
  if (crop_model == 'wheat1') {
    data("GDD_wheat1")
    base_temperature = 0
    max_temperature = 35
    return(list(GDD_wheat1,base_temperature,max_temperature))
  }
  if (crop_model == 'wheat2') {
    data("GDD_wheat2")
    base_temperature = 0
    max_temperature = 35
    return(list(GDD_wheat2,base_temperature,max_temperature))
  }
  if (crop_model == 'hardwheatUS') {
    data("GDD_hardredwheatUS")
    base_temperature = 0
    max_temperature = 35
    return(list(GDD_hardredwheatUS,base_temperature,max_temperature))
  }
  if (crop_model == 'barley') {
    data("GDD_barley")
    base_temperature = 0
    max_temperature = 35
    return(list(GDD_barley,base_temperature,max_temperature))
  }
}