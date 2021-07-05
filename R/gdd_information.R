gdd_information <- function(crop_model){
  if (crop_model == 'maizehybrid2700'){
    data("GDD_maize2700")
    base_temperature = 10
    return(list(GDD_maize2700,base_temperature))
  }
}