devtools::load_all()

data("geno_indica")
data("map_indica")
data("pheno_indica")
data("info_environments_indica")
data("climate_variables_indica")

METdata_indica <-
  create_METData(
    geno = geno_indica,
    pheno = pheno_indica,
    climate_variables = climate_variables_indica,
    compute_climatic_ECs = F,
    info_environments = info_environments_indica,
    map = map_indica,
    path_to_save = '~/scripts/INDICA'
  )


# Trait GY

start <- Sys.time()
print(start)


rescv0_2 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'GY',
  prediction_method = 'stacking_reg_1',
  use_selected_markers = F,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  kernel_G = 'linear',
  include_env_predictors = T,
  save_processing  = T,
  seed = 100,
  path_folder = '~/scripts/INDICA/stacking_reg_1/cv0'
)


end_time <- Sys.time()
print(end_time - start)
saveRDS(rescv0_2,
        '~/scripts/INDICA/stacking_reg_1/cv0/rescv0_gy.RDS')


# Trait PHR

start <- Sys.time()
print(start)


rescv0_3 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'PHR',
  prediction_method = 'stacking_reg_1',
  use_selected_markers = F,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  kernel_G = 'linear',
  include_env_predictors = T,
  save_processing  = T,
  seed = 100,
  path_folder = '~/scripts/INDICA/stacking_reg_1/cv0'
)


end_time <- Sys.time()
print(end_time - start)
saveRDS(rescv0_3,
        '~/scripts/INDICA/stacking_reg_1/cv0/rescv0_phr.RDS')


# Trait GC

start <- Sys.time()
print(start)


rescv0_4 <- predict_trait_MET_cv(
  METData = METdata_indica,
  trait = 'GC',
  prediction_method = 'stacking_reg_1',
  use_selected_markers = F,
  lat_lon_included = F,
  year_included = F,
  cv_type = 'cv0',
  cv0_type = 'leave-one-year-out',
  kernel_G = 'linear',
  include_env_predictors = T,
  save_processing  = T,
  seed = 100,
  path_folder = '~/scripts/INDICA/stacking_reg_1/cv0'
)


end_time <- Sys.time()
print(end_time - start)
saveRDS(rescv0_4,
        '~/scripts/INDICA/stacking_reg_1/cv0/rescv0_gc.RDS')
