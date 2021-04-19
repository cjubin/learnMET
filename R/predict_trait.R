#' Trait prediction based on genetic and environmental data
#' This function assumes that METData$pheno contains all
#'
#' @param METData. An object created by the initial function of the package,
#' "create_METData.R"
#'
#' @param trait. \code{character} Name of the trait under study for which a
#'
#' @param
#'
#'
#' @return
#'
#'


predict_trait_MET_cv <- function(METData,
                                 trait,
                                 method = 'xgboost',
                                 use_selected_markers = T,
                                 lat_lon_included = T,
                                 year_included = F,
                                 cv_type = c('cv0', 'cv1', 'cv2'),
                                 cv0_type = c('leave-one-environment-out',
                                              'leave-one-site-out',
                                              'leave-one-year-out'),
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 list_env_predictors = colnames(METData$env_data)[colnames(METData$env_data) %notin%
                                                                                    c('IDenv', 'year', 'location', 'longitude', 'latitude')]
                                 ,...) {
 
  geno = METData$geno
  pheno = METData$pheno
  env_predictors = METData$env_data
  
  # Select the genotypic data to use if a subset of selected markers was created
  # in a previous step.
  
  if (use_selected_markers == T) {
    geno = geno[, colnames(METData$geno) %in% METData$selected_markers]
  }
  
  # Merge in same data pheno and geno data, with environmental covariates
  
  geno$geno_ID = row.names(geno)
  pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)
  pheno <- merge(pheno, env_predictors, by = 'IDenv', all.x = T)
  
  if (lat_lon_included) {
    list_env_predictors = c(list_env_predictors, 'latitude', 'longitude')
  }
  
  if (year_included) {
    list_env_predictors = c(list_env_predictors, 'year')
  }
  
   # Select trait from the phenotypic file
  
  list_predictors_geno <- colnames(geno)
  
  pheno <-
    pheno[, c(trait, list_predictors_geno, list_env_predictors)]
  
  
  # Create cross-validation random splits according to the type of selected CV
  
  if (cv_type == 'cv1') {
    splits <-
      predict_cv1(pheno_all_data = pheno,
                  nb_folds = nb_folds_cv1,
                  reps = repeats_cv1)
  }
  
  if (cv_type == 'cv2') {
    splits <-
      predict_cv2(pheno_all_data = pheno,
                  nb_folds = nb_folds_cv2,
                  reps = repeats_cv2)
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #Create the cross-validation random splits for hyperparameter optimization
  
  cv_splits <-
    rsample::vfold_cv(pheno, folds = nb_folds_cv, repeats = reps)
  
  # Define the type of model to tune
  
  mod <- parsnip::linear_reg(penalty = tune(),
                             mixture = tune()) %>%
    parsnip::set_engine("glmnet")
  
  
  # Define the predictors and outcome variables.
  # Remove predictors with null varaiance.
  # Center and scale the variables (standardization)
  
  rec <- recipes::recipe( ~ ., data = pheno) %>%
    recipes::update_role(all_of(trait), new_role = 'outcome') %>%
    recipes::update_role(all_of(list_predictors), new_role = 'predictor') %>%
    recipes::step_nzv(all_numeric()) %>%
    recipes::step_normalize(all_numeric())
  
  
}


#@reps:number of k - folds random partitions

