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
    return(list(GDD_maize1700,base_temperature))
  }
  if (crop_model == 'hardwheatUS') {
    data("GDD_hardredwheatUS")
    base_temperature = 0
    return(list(GDD_hardredwheatUS,base_temperature))
    
  }
}