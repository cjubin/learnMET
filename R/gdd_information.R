gdd_information <- function(crop_model){
  if (crop_model == 'maizehybrid1700'){
    data("GDD_maize1700")
    base_temperature = 10
    return(list(GDD_maize1700,base_temperature))
  }
}