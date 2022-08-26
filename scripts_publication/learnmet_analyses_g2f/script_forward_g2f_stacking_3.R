setwd('~/Data/PackageMLpredictions/learnMET')
devtools::load_all()
library(dplyr)

data(geno_G2F)
#geno_G2F <- geno_G2F[,1:4500]
data(pheno_G2F)
data(map_G2F)
#map_G2F <- map_G2F[map_G2F$marker%in%colnames(geno_G2F),]
data(info_environments_G2F)
data(soil_G2F)


METdata_G2F_training <-
  create_METData(
    geno = geno_G2F,
    pheno = pheno_G2F[pheno_G2F$year%in%c(2014,2015,2016),],
    map = map_G2F,
    climate_variables = NULL,
    compute_climatic_ECs = TRUE,
    et0 =T,
    info_environments = info_environments_G2F[info_environments_G2F$year%in%c(2014,2015,2016),],
    soil_variables = soil_G2F[soil_G2F$year%in%c(2014,2015,2016),],
    path_to_save = "~/results_benchmarking_g2f/results_g2f_forward_3"
  )


METdata_G2F_new <-
  create_METData(
    geno = geno_G2F,
    pheno = as.data.frame(pheno_G2F[pheno_G2F$year%in%2017,] %>% dplyr::select(-pltht,-yld_bu_ac,-earht)),
    map = map_G2F,
    et0=T,
    climate_variables = NULL,
    compute_climatic_ECs = TRUE,
    info_environments = info_environments_G2F[info_environments_G2F$year%in%2017,],
    soil_variables = soil_G2F[soil_G2F$year%in%2017,],
    path_to_save = "~/results_benchmarking_g2f/results_g2f_forward_3",
    as_test_set = T
  )


start <- Sys.time()
print(start)
rescv0_2 <- predict_trait_MET(
  METData_training = METdata_G2F_training,
  METData_new = METdata_G2F_new,
  trait = 'yld_bu_ac',
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
  path_folder = '~/results_benchmarking_g2f/results_g2f_forward_3/stacking_reg_3'
)
end_time <- Sys.time()
print(end_time - start)







pheno_G2F$IDenv <- paste0(pheno_G2F$location,'_',pheno_G2F$year)


met_pred <- readRDS("~/results_benchmarking_g2f/results_g2f_forward_3/stacking_reg_3/met_pred.RDS")
pred_df <- met_pred$list_results[[1]]$predictions_df
pred_df <- plyr::join(pred_df %>% dplyr::select(-yld_bu_ac,-year,-location),  as.data.frame(pheno_G2F[pheno_G2F$year%in%2017,]),by=c('IDenv','geno_ID') )

pred_df %>% group_by(IDenv) %>% summarise(cor = cor(yld_bu_ac, .pred))
pred_df %>% group_by(IDenv) %>% yardstick::rmse(yld_bu_ac, .pred)

