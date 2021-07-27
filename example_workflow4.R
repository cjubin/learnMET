setwd("/home/uni08/jubin1/Data/PackageMLpredictions/learnMET")

devtools::load_all()
library(ggplot2)
library(assertive.datetimes)
library(purrr)
library(ggrepel)
library(nasapower)
library(tidymodels)
library(keras)
library(vip)
library(purrr)
library(tidytext)
METdata_indica <-
  create_METData(
    geno = geno_indica,
    pheno = pheno_indica,
    env_data = env_data_indica,
    compute_ECs = F,
    info_environments = info_environments_indica,
    map = map_indica
  )

METdata_indica$geno <- METdata_indica$geno[, 1:5000]

start_time <- Sys.time()
rescv0_1 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'PH',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = F,
  geno_information = 'SNPs',
  num_pcs = 300,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  kernel_G = 'linear',
  include_env_predictors = T,
  list_env_predictors = NULL,
  save_processing  = T,
  seed = 100,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/cv0'
)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(rescv0_1,'/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/rescv0_1.RDS')




start_time <- Sys.time()
rescv0_1 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'PHR',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = F,
  geno_information = 'SNPs',
  num_pcs = 300,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  kernel_G = 'linear',
  include_env_predictors = T,
  list_env_predictors = NULL,
  save_processing  = T,
  seed = 100,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/cv0'
)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(rescv0_1,'/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/rescv0_1.RDS')


start_time <- Sys.time()
rescv0_1 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'GY',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = F,
  geno_information = 'SNPs',
  num_pcs = 300,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  kernel_G = 'linear',
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  save_processing  = T,
  seed = 100,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/cv0'
)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(rescv0_1,'/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/rescv0_1.RDS')



start_time <- Sys.time()
rescv0_1 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'GC',
  method_processing = 'svm_stacking_reg',
  use_selected_markers = F,
  geno_information = 'SNPs',
  num_pcs = 300,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  nb_folds_cv1 = 3,
  repeats_cv1 = 2,
  kernel_G = 'linear',
  nb_folds_cv2 = 5,
  repeats_cv2 = 50,
  include_env_predictors = T,
  list_env_predictors = NULL,
  save_processing  = T,
  seed = 100,
  path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/cv0'
)
end_time <- Sys.time()
print(end_time - start_time)
saveRDS(rescv0_1,'/home/uni08/jubin1/Data/PackageMLpredictions/learnMET/INDICA/svm_noGE/rescv0_1.RDS')
