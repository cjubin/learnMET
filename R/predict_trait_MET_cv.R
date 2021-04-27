#' Trait prediction based on SNP and environmental data
#' This function should be used to assess the predictive ability according to
#' various cross-validation schemes.
#'
#' @param METData. An object created by the initial function of the package,
#' "create_METData.R"
#'
#' @param trait. \code{character} Name of the trait under study for which a
#' 
#' @param use_selected_markers \code{Logical} Whether to use a subset of markers
#' obtained from a previous step (see function select_markers()).
#'
#' 
#'
#'
#' 
#'
#'


predict_trait_MET_cv <- function(METData,
                                 trait,
                                 method = 'xgboost',
                                 use_selected_markers = T,
                                 geno_information = c('SNPs','PCS','PCs_G'),
                                 num_pcs = 200,
                                 lat_lon_included = T,
                                 year_included = F,
                                 cv_type = c('cv0', 'cv1', 'cv2'),
                                 cv0_type = c('leave-one-environment-out',
                                              'leave-one-site-out',
                                              'leave-one-year-out',
                                              'forward-prediction'),
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 list_env_predictors = colnames(METData$env_data)[colnames(METData$env_data) %notin%
                                                                                    c('IDenv', 'year', 'location', 'longitude', 'latitude')]
                                 ,...) {
  
  geno = METData$geno
  env_predictors = METData$env_data
  
  # Select phenotypic data for the trait under study and remove NA in phenotypes
  
  pheno = METData$pheno[,c("geno_ID","year" ,"location","IDenv",trait)][complete.cases(METData$pheno[,c("geno_ID","year" ,"location","IDenv",trait)]),]
  
  
  # Create cross-validation random splits according to the type of selected CV
  
  if (cv_type == 'cv1') {
    splits <-
      predict_cv1(pheno_data = pheno,
                  nb_folds = nb_folds_cv1,
                  reps = repeats_cv1)
  }
  
  if (cv_type == 'cv2') {
    splits <-
      predict_cv2(pheno_data = pheno,
                  nb_folds = nb_folds_cv2,
                  reps = repeats_cv2)
  }
  
  
  
  if (cv_type == 'cv0') {
    splits <-
      predict_cv0(pheno_data = pheno,
                  type = cv0_type)
  }
  
  ###############################
  ###############################
  
  
  ## USE OF GENOTYPIC DATA ##
  
  rec <- recipe(yld_bu_ac ~ . ,
                data = training) %>%
    update_role(yld_bu_ac, new_role = 'outcome') %>%
    update_role(year, new_role = "id variable") %>%
    update_role(Year_Exp, new_role = "id variable") %>%
    update_role(-yld_bu_ac, -year,-Year_Exp, new_role = 'predictor') %>%
    step_nzv(all_predictors(),-starts_with('PC')) %>%
    step_normalize(all_numeric(), -all_outcomes(),-starts_with('PC'))
  
  
  if (geno_information == 'PCs'){
    
    # Processing of PCs: apply transformations calculated on the training set, 
    # on test set.
    
    apply_pca <-function(split){
      
      geno$geno_ID = row.names(geno)
      
      col_to_keep = colnames(geno)
      
      geno_training = merge(split[[1]], geno, by = 'geno_ID', all.x = T)[,col_to_keep]
      geno_training = unique(geno_training)
      
      geno_test =  merge(split[[2]], geno, by = 'geno_ID', all.x = T)[,col_to_keep]
      geno_test = unique(geno_test)
      
      
        rec <- recipes::recipe(geno_ID ~ . ,
                               data = geno_training) %>%
          update_role(geno_ID, new_role = 'outcome') %>%
          recipes::step_pca(all_predictors(), num_comp = num_pcs,options = list(center=T,scale.=T))
        
        norm_obj <- prep(rec, training = geno_training)
        
        training_pca <- juice(norm_obj)
      
         
      test_pca <- bake(norm_obj, geno_test)
      test_pca$geno_ID=geno_test$geno_ID
      
      pc_values<-rbind(training_pca,test_pca)
      pc_values<-unique(pc_values)
      
    
  }
  
  # Select the SNPs to use if a subset of selected markers was created
  # in a previous step.
  
  if (use_selected_markers == T) {
    SNPs = geno[, colnames(METData$geno) %in% METData$selected_markers]
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
  
  # Process the training and test sets
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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

