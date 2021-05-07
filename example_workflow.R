setwd("/home/uni08/jubin1/Data/PackageMLpredictions/multienvtML")

devtools::load_all()
data(geno_indica)
data(map_indica)
data(pheno_indica)
data(info_environments_indica)
data(env_data_indica)

#library(purrr)
METdata_indica <-
  create_METData(
    geno = geno_indica,
    pheno = pheno_indica,
    env_data = env_data_indica,
    unique_EC_by_geno = F,
    compute_ECs = F,
    info_environments = info_environments_indica,
    map = map_indica
  )

METdata_indica$geno<-METdata_indica$geno[,1:5000]

METdata_indica2 <- select_markers(
  METData = METdata_indica,
  trait = 'PH',
  method_marker_effects = 'FarmCPU',
  method_selection_EN  = c('only_variance_across_env'),
  size_subset_most_variable_markers = 200,
  size_top_markers_by_env = 400,
  plot_penalty_regression_coefficients = F,
  plot_gwas = T,
  path_save_plot =  "/home/uni08/jubin1/Data/PackageMLpredictions/plots",
  path_save_results =  "/home/uni08/jubin1/Data/PackageMLpredictions/try_indica"
)

#saveRDS(METdata_indica2,'/home/uni08/jubin1/Data/PackageMLpredictions/try_indica/METdata_indica2farmCPU.RDS')

res0 <- predict_trait_MET_cv(METdata_indica2,
                            trait='PH',
                            method = 'xgboost',
                            use_selected_markers = T,
                            geno_information = c('PCs'),
                            num_pcs = 200,
                            lat_lon_included = T,
                            year_included = T,
                            cv_type = c('cv1'),
                            cv0_type = c(
                              'leave-one-year-out'
                            ),
                            nb_folds_cv1 = 3,
                            repeats_cv1 = 2,
                            nb_folds_cv2 = 5,
                            repeats_cv2 = 50,
                            include_env_predictors = T,
                            list_env_predictors = NULL,
                            plot_PA = T,
                            path_plot_PA = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/'
)


res <- predict_trait_MET_cv(METdata_indica2,
                                 trait='PH',
                                 method = 'xgboost',
                                 use_selected_markers = T,
                                 geno_information = c('PCs'),
                                 num_pcs = 200,
                                 lat_lon_included = T,
                                 year_included = T,
                                 cv_type = c('cv0'),
                                 cv0_type = c(
                                   'leave-one-year-out'
                                 ),
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 include_env_predictors = T,
                                 list_env_predictors = NULL,
                                 plot_PA = T,
                                 path_plot_PA = '/home/uni08/jubin1/Data/PackageMLpredictions/plots/'
                                 )

saveRDS(res,'/home/uni08/jubin1/Data/PackageMLpredictions/try_indica/predictions_cv0.RDS')



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