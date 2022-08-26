setwd('~/Data/PackageMLpredictions/learnMET')
devtools::load_all()
library(dplyr)

data("geno_indica")
data("map_indica")
data("pheno_indica")
data("info_environments_indica")
data("climate_variables_indica")


METdata_indica_training <-
  create_METData(
    geno = geno_indica,
    pheno = pheno_indica[pheno_indica$year%in%c(2010,2011),],
    climate_variables = climate_variables_indica[climate_variables_indica$year%in%c(2010,2011),],
    compute_climatic_ECs = F,
    map = map_indica,
    et0 =T,
    info_environments = info_environments_indica[info_environments_indica$year%in%c(2010,2011),],
    path_to_save = "~/results_benchmarking_rice_indicaPHR",
  )


METdata_indica_new <-
  create_METData(
    geno = geno_indica,
    pheno = as.data.frame(pheno_indica[pheno_indica$year%in%2012,] %>% dplyr::select(-PHR,-GY,-/GC)),
    map = map_indica,
    et0=T,
    climate_variables = climate_variables_indica[climate_variables_indica$year%in%2012,],
    compute_climatic_ECs = F,
    info_environments = info_environments_indica[info_environments_indica$year%in%2012,],
    path_to_save = "~/results_benchmarking_rice_indicaPHR",
    as_test_set = T
  )


start <- Sys.time()
print(start)
rescv0_2 <- predict_trait_MET(
  METData_training = METdata_indica_training,
  METData_new = METdata_indica_new,
  trait = 'PHR',
  prediction_method = 'xgb_reg_1',
  use_selected_markers = F,
  save_model = T,
  lat_lon_included = F,
  year_included = F,
  num_pcs = 100,
  include_env_predictors = T,
  save_splits = T,
  seed = 100,
  save_processing = T,
  path_folder = '~/results_benchmarking_rice_indicaPHR/xgbreg/fw2'
)
end_time <- Sys.time()
print(end_time - start)



pheno_indica$IDenv <- paste0(pheno_indica$location,'_',pheno_indica$year)


pred_df <- rescv0_2$list_results[[1]]$predictions_df
pred_df <- plyr::join(pred_df %>% dplyr::select(-GY,-/GC,-PHR,-year,-location),  as.data.frame(pheno_indica[pheno_indica$year%in%2012,]),by=c('IDenv','geno_ID') )

pred_df %>% group_by(IDenv) %>% summarise(cor = cor(PHR, .pred))






















