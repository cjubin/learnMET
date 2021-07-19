setwd("/home/uni08/jubin1/Data/PackageMLpredictions/multienvtML")

devtools::load_all()
library(ggplot2)
library(assertive.datetimes)
library(purrr)
library(ggrepel)
library(nasapower)
library(tidymodels)
library(keras)
library(vip)
data(geno_G2F)
data(pheno_G2F)
data(map_G2F)
data(info_environments_G2F)

pheno_G2F_1 <-
  pheno_G2F[pheno_G2F$location %in% c('Georgetown', 'Columbia'), ]
geno_G2F_1 <- geno_G2F[row.names(geno_G2F) %in% pheno_G2F_1$geno_ID, ]

info_environments_G2F_1 <-
  info_environments_G2F[info_environments_G2F$location %in% c('Georgetown', 'Columbia'), ]
METdata_g2f <-
  create_METData(
    geno = geno_G2F_1,
    pheno = pheno_G2F_1,
    map = map_G2F,
    env_data = NULL,
    compute_climatic_ECs = TRUE,
    info_environments = info_environments_G2F_1,
    crop_model = 'maizehybrid1700'
  )


METdata_g2f$geno <- METdata_g2f$geno[, 1:5000]

pheno_new <- pheno_G2F[pheno_G2F$location %in% c('CollegeStation')&pheno_G2F$year==2015, ]
pheno_new <- pheno_new %>% dplyr::select(geno_ID, year, location)
unique_ID <- unique(pheno_new$geno_ID)
pheno_2014 <- cbind(geno_ID=unique_ID,year='2014',location='CollegeStation') 
pheno_2013 <- cbind(geno_ID=unique_ID,year='2013',location='CollegeStation')
pheno_2016 <- cbind(geno_ID=unique_ID,year='2016',location='CollegeStation')
pheno_2017 <- cbind(geno_ID=unique_ID,year='2017',location='CollegeStation')
pheno_2018 <- cbind(geno_ID=unique_ID,year='2018',location='CollegeStation')
pheno_2019 <- cbind(geno_ID=unique_ID,year='2019',location='CollegeStation')
pheno_2020 <- cbind(geno_ID=unique_ID,year='2020',location='CollegeStation')
pheno_aurora0 <- cbind(geno_ID=unique_ID,year='2018',location='Aurora')
pheno_aurora1 <- cbind(geno_ID=unique_ID,year='2019',location='Aurora')
pheno_aurora2 <- cbind(geno_ID=unique_ID,year='2020',location='Aurora')
pheno_new <- rbind(pheno_2013,pheno_2014,pheno_new,pheno_2016,pheno_2017,pheno_2018,pheno_2019,pheno_2020,pheno_aurora0,pheno_aurora1,pheno_aurora2)

geno_new <- geno_G2F[which(row.names(geno_G2F)%in%pheno_new$geno_ID), ]
info_environments_to_predict <-
  info_environments_G2F[info_environments_G2F$location %in% c('CollegeStation'), ]
info_environments_to_predict <- rbind(info_environments_to_predict, c(2013,'CollegeStation',unique(info_environments_to_predict$longitude)[1],unique(info_environments_to_predict$latitude)[1], stringr::str_replace(info_environments_to_predict$planting.date[1],"2014",'2013'), stringr::str_replace(info_environments_to_predict$harvest.date[1],"2014",'2013')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2017,'CollegeStation',unique(info_environments_to_predict$longitude)[1],unique(info_environments_to_predict$latitude)[1], stringr::str_replace(info_environments_to_predict$planting.date[1],"2014",'2017'), stringr::str_replace(info_environments_to_predict$harvest.date[1],"2014",'2017')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2018,'CollegeStation',unique(info_environments_to_predict$longitude)[1],unique(info_environments_to_predict$latitude)[1], stringr::str_replace(info_environments_to_predict$planting.date[1],"2014",'2018'), stringr::str_replace(info_environments_to_predict$harvest.date[1],"2014",'2018')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2019,'CollegeStation',unique(info_environments_to_predict$longitude)[1],unique(info_environments_to_predict$latitude)[1], stringr::str_replace(info_environments_to_predict$planting.date[1],"2014",'2019'), stringr::str_replace(info_environments_to_predict$harvest.date[1],"2014",'2019')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2020,'CollegeStation',unique(info_environments_to_predict$longitude)[1],unique(info_environments_to_predict$latitude)[1], stringr::str_replace(info_environments_to_predict$planting.date[1],"2014",'2020'), stringr::str_replace(info_environments_to_predict$harvest.date[1],"2014",'2020')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2018,'Aurora',-76.65000,42.73000, stringr::str_replace(info_environments_to_predict$planting.date[1],as.character(info_environments_to_predict$planting.date[1]),'2018-05-15'), stringr::str_replace(info_environments_to_predict$harvest.date[1],as.character(info_environments_to_predict$harvest.date[1]),'2018-11-25')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2019,'Aurora',-76.65000,42.73000, stringr::str_replace(info_environments_to_predict$planting.date[1],as.character(info_environments_to_predict$planting.date[1]),'2019-05-15'), stringr::str_replace(info_environments_to_predict$harvest.date[1],as.character(info_environments_to_predict$harvest.date[1]),'2019-11-25')))
info_environments_to_predict <- rbind(info_environments_to_predict, c(2020,'Aurora',-76.65000,42.73000, stringr::str_replace(info_environments_to_predict$planting.date[1],as.character(info_environments_to_predict$planting.date[1]),'2020-05-15'), stringr::str_replace(info_environments_to_predict$harvest.date[1],as.character(info_environments_to_predict$harvest.date[1]),'2020-11-25')))
class(pheno_new$year)<-'numeric'
class(info_environments_to_predict$year)<-'numeric'
info_environments_to_predict$longitude<-as.numeric(as.character(info_environments_to_predict$longitude))
info_environments_to_predict$latitude<-as.numeric(as.character(info_environments_to_predict$latitude))


METdata_to_predict_CollegeStation <-
  add_new_METData(
    geno_new = geno_new,
    METData_training = METdata_g2f,
    pheno_new = pheno_new,
    compute_climatic_ECs = TRUE,
    info_environments_to_predict = info_environments_to_predict,
    method_ECs_intervals = 'fixed_nb_windows_across_env'
  )


predicted_new_data_with_singlemarkers <-
  predict_trait_MET(
    METData_training = METdata_g2f,
    METData_new = METdata_to_predict_CollegeStation,
    method_processing = 'xgb_reg',
    use_selected_markers = T,
    geno_information = 'PCs',
    num_pcs = 50,
    trait = 'pltht',
    lat_lon_included = F,
    path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/to_predict',
    vip_plot = F
  )









































predicted_new_data_CollegeStation <-
  predict_trait_MET(
    METData_training = METdata_g2f,
    METData_new = METdata_to_predict_CollegeStation,
    method_processing = 'xgb_reg',
    geno_information = 'PCs',
    num_pcs = 50,
    trait = 'pltht',
    lat_lon_included = F,
    path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/to_predict',
    vip_plot = F
  )






pheno_new <- pheno_G2F[pheno_G2F$location %in% c('CollegeStation'), ]
pheno_new <- pheno_new %>% select(geno_ID, year, location)
geno_new <- geno_G2F[row.names(geno_G2F) %in% pheno_new$geno_ID, ]
info_environments_to_predict <-
  info_environments_G2F[info_environments_G2F$location %in% c('CollegeStation'), ]

METdata_to_predict <-
  add_new_METData(
    geno_new = geno_new,
    METData_training = METdata_g2f,
    pheno_new = pheno_new,
    compute_climatic_ECs = TRUE,
    info_environments_to_predict = info_environments_to_predict,
    crop_model = 'maizehybrid1700'
  )

predicted_new_data <-
  predict_trait_MET(
    METData_training = METdata_g2f,
    METData_new = METdata_to_predict,
    list_env_predictors = colnames(METdata_to_predict$env_data)[4:12],
    method_processing = 'xgb_reg',
    geno_information = 'PCs',
    num_pcs = 50,
    trait = 'pltht',
    lat_lon_included = F,
    path_folder = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/g2f/pltht/to_predict',
    vip_plot = F
  )





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
  vip_plot = F
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
  vip_plot = F
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
