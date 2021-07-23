gdd_information <- function(crop_model){
  if (crop_model == 'maizehybrid1700'){
    data("GDD_maize1700")
    base_temperature = 10
    return(list(GDD_maize1700,base_temperature))
  }
  if (crop_model == 'hardwheatUS'){
    data("GDD_hardredwheatUS")
    base_temperature = 0
    return(list(GDD_hardredwheatUS,base_temperature))
    
  }
}