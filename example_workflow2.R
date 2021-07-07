setwd("/home/uni08/jubin1/Data/PackageMLpredictions/multienvtML")

devtools::load_all()
library(ggplot2)
library(assertive.datetimes)
library(purrr)
library(ggrepel)
data(geno_G2F)
data(pheno_G2F)
data(map_G2F)
data(info_environments_G2F)
pheno_G2F <- pheno_G2F[pheno_G2F$location%in%c('Georgetown','Columbia'),]
info_environments_G2F <- info_environments_G2F[info_environments_G2F$location%in%c('Georgetown','Columbia'),]
METdata_g2f <- create_METData(geno=geno_G2F,pheno=pheno_G2F,map=map_G2F,env_data = NULL,compute_ECs = TRUE,info_environments = info_environments_G2F,crop_model='maizehybrid1700')


METdata_g2f$geno <- METdata_g2f$geno[, 1:10000]


#saveRDS(METdata_indica2,'/home/uni08/jubin1/Data/PackageMLpredictions/try_indica/METdata_indica2farmCPU.RDS')
rescv0_g2f <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 100,
  lat_lon_included = T,
  year_included = F,
  cv_type = c('cv0'),
  cv0_type = c('leave-one-year-out'),
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  plot_PA = T,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/cv0'
)

rescv1 <- predict_trait_MET_cv(
  METData = METdata_indica2,
  trait = 'PH',
  method_processing = 'xgb_reg',
  use_selected_markers = T,
  geno_information = c('PCs'),
  num_pcs = 200,
  lat_lon_included = T,
  year_included = F,
  cv_type = c('cv1'),
  cv0_type = c('leave-one-year-out'),
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  plot_PA = T,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/cv1'
)

res_cv2 <- predict_trait_MET_cv(
  METData = METdata_indica2,
  trait = 'PH',
  method_processing = 'xgb_reg',
  use_selected_markers = T,
  geno_information = c('PCs'),
  num_pcs = 200,
  lat_lon_included = T,
  year_included = T,
  cv_type = 'cv2',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 3,
  repeats_cv2 = 2,
  include_env_predictors = T,
  list_env_predictors = NULL,
  plot_PA = T,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/cv2'
)

res_svm <- predict_trait_MET_cv(
  METData = METdata_indica2,
  trait = 'PH',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = T,
  geno_information = c('PCs'),
  num_pcs = 200,
  lat_lon_included = T,
  year_included = T,
  cv_type = 'cv2',
  cv0_type = c('leave-one-environment-out'),
  nb_folds_cv1 = 5,
  repeats_cv1 = 50,
  nb_folds_cv2 = 3,
  repeats_cv2 = 2,
  include_env_predictors = T,
  list_env_predictors = NULL,
  plot_PA = T,
  filename_plot_PA = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/',
  save_processing = F,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/',
  kernel_G = 'rbf',
  kernel_GE = 'rbf',
  kernel_E = 'rbf',
  variable_importance = T
)

saveRDS(
  res,
  '/home/uni08/jubin1/Data/PackageMLpredictions/try_indica/predictions_cv0.RDS'
)



## OR
METdata_indica3 <- select_markers(
  METdata_indica,
  trait = 'PH',
  method_marker_effects = 'elasticnet',
  method_selection = c('effect_size_per_env'),
  size_subset_most_variable_markers = 200,
  size_top_markers_by_env = 400,
  plot_penalty_regression_coefficients = F,
  plot_gwas = T,
  path_save_plot =  "/home/uni08/jubin1/Data/PackageMLpredictions/plots",
  path_save_results =  "/home/uni08/jubin1/Data/PackageMLpredictions/try_indica"
)