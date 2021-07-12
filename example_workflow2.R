setwd("/home/uni08/jubin1/Data/PackageMLpredictions/multienvtML")

devtools::load_all()
library(ggplot2)
library(assertive.datetimes)
library(purrr)
library(ggrepel)f
library(nasapower)
library(tidymodels)
library(keras)
data(geno_G2F)
data(pheno_G2F)
data(map_G2F)
data(info_environments_G2F)
pheno_G2F <- pheno_G2F[pheno_G2F$location%in%c('Georgetown','Columbia'),]
info_environments_G2F <- info_environments_G2F[info_environments_G2F$location%in%c('Georgetown','Columbia'),]
METdata_g2f <- create_METData(geno=geno_G2F,pheno=pheno_G2F,map=map_G2F,env_data = NULL,compute_ECs = TRUE,info_environments = info_environments_G2F,crop_model='maizehybrid1700')


METdata_g2f$geno <- METdata_g2f$geno[, 1:12000]


#saveRDS(METdata_indica2,'/home/uni08/jubin1/Data/PackageMLpredictions/try_indica/METdata_indica2farmCPU.RDS')
rescv0_g2f_pltht_0 <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/cv0'
)

rescv0_g2f_pltht_1 <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-site-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/cv0'
)

rescv0_g2f_pltht_2 <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = c('forward-prediction'),
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/cv0'
)


rescv2_g2f_pltht <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv2',
  cv0_type = c('forward-prediction'),
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 4,
  repeats_cv2 = 2,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/cv2'
)

rescv1_g2f <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'xgb_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv1',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 4,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/cv1'
)


METdata_g2f_markers <- select_markers(
  METData = METdata_indica,
  trait = 'pltht',
  method_marker_effects = 'FarmCPU',
  method_selection_EN  = c('only_variance_across_env'),
  size_subset_most_variable_markers = 200,
  size_top_markers_by_env = 400,
  plot_penalty_regression_coefficients = F,
  plot_gwas = T,
  path_save_plot =  "/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f",
  path_save_results =  "/home/uni08/jubin1/Data/PackageMLpredictions/try_g2f"
)

res_cv0_g2f0_svm <- predict_trait_MET_cv(
  METData = METdata_g2f_markers,
  trait = 'pltht',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = T,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 3,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/svm/pcs/cv0',
  vip_plot =F
)

res_cv0_g2f0_svm_snps <- predict_trait_MET_cv(
  METData = METdata_g2f_markers,
  trait = 'pltht',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = T,
  geno_information = c('SNPs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 3,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/svm/snps/cv0',
  vip_plot =F
)

rescv0_g2f_pltht_0 <- predict_trait_MET_cv(
  METData = METdata_g2f,
  trait = 'pltht',
  method_processing = 'DL_reg',
  use_selected_markers = F,
  geno_information = c('PCs'),
  num_pcs = 150,
  lat_lon_included = T,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/DL/cv0'
)


