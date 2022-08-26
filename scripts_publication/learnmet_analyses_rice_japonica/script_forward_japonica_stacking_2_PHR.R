setwd('~/Data/PackageMLpredictions/learnMET')
devtools::load_all()
library(dplyr)

data("geno_japonica")
data("map_japonica")
data("pheno_japonica")
data("info_environments_japonica")
data("climate_variables_japonica")


METdata_japonica_training <-
  create_METData(
    geno = geno_japonica,
    pheno = pheno_japonica[pheno_japonica$year%in%c(2009,2010),],
    climate_variables = climate_variables_japonica[climate_variables_japonica$year%in%c(2009,2010),],
    compute_climatic_ECs = F,
    map = map_japonica,
    et0 =T,
    info_environments = info_environments_japonica[info_environments_japonica$year%in%c(2009,2010),],
    path_to_save = "~/results_benchmarking_rice_japonica/PHR",
  )


METdata_japonica_new <-
  create_METData(
    geno = geno_japonica,
    pheno = as.data.frame(pheno_japonica[pheno_japonica$year%in%2011,] %>% dplyr::select(-PHR,-GY,-GC)),
    map = map_japonica,
    et0=T,
    climate_variables = climate_variables_japonica[climate_variables_japonica$year%in%2011,],
    compute_climatic_ECs = F,
    info_environments = info_environments_japonica[info_environments_japonica$year%in%2011,],
    path_to_save = "~/results_benchmarking_rice_japonica/PHR",
    as_test_set = T
  )


start <- Sys.time()
print(start)
rescv0_2 <- predict_trait_MET(
  METData_training = METdata_japonica_training,
  METData_new = METdata_japonica_new,
  trait = 'PHR',
  prediction_method = 'stacking_reg_3',
  use_selected_markers = F,
  save_model = F,
  lat_lon_included = F,
  year_included = F,
  num_pcs = 100,
  include_env_predictors = T,
  save_splits = T,
  seed = 100,
  save_processing = T,
  path_folder = '~/results_benchmarking_rice_japonica/PHR/stacking/fw2'
)
end_time <- Sys.time()
print(end_time - start)










pheno_japonica$IDenv <- paste0(pheno_japonica$location,'_',pheno_japonica$year)


pred_df <- rescv0_2$list_results[[1]]$predictions_df
pred_df <- plyr::join(pred_df %>% dplyr::select(-GY,-GC,-PHR,-year,-location),  as.data.frame(pheno_japonica[pheno_japonica$year%in%2011,]),by=c('IDenv','geno_ID') )

pred_df %>% group_by(IDenv) %>% summarise(cor = cor(PHR, .pred))
















