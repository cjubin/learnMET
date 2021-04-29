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
                                 geno_information = c('SNPs', 'PCs', 'PCs_G'),
                                 num_pcs = 200,
                                 lat_lon_included = T,
                                 year_included = F,
                                 cv_type = c('cv0', 'cv1', 'cv2'),
                                 cv0_type = c(
                                   'leave-one-environment-out',
                                   'leave-one-site-out',
                                   'leave-one-year-out',
                                   'forward-prediction'
                                 ),
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 include_env_predictors = T,
                                 list_env_predictors = NULL
                                 ,
                                 ...) {
  geno = METData$geno
  
  # If no specific list of environmental predictors is provided by the user:
  # all of the environmental predictors present in METData$env_data are used as
  # predictors.
  
  if (is.null(list_env_predictors) & include_env_predictors) {
    list_env_predictors = colnames(METData$env_data)[colnames(METData$env_data) %notin%
                                                       c('IDenv', 'year', 'location', 'longitude', 'latitude')]
    
    
  }
  
  env_predictors = METData$env_data
  
  
  # Select phenotypic data for the trait under study and remove NA in phenotypes
  
  pheno = METData$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)][complete.cases(METData$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)]), ]
  
  
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
  
  if (geno_information == 'PCs') {
    # Processing of PCs: apply transformations calculated on the training set,
    # on test set.
    
    apply_pca <- function(split) {
      geno$geno_ID = row.names(geno)
      
      col_to_keep = colnames(geno)
      
      geno_training = merge(split[[1]], geno, by = 'geno_ID', all.x = T)[, col_to_keep]
      geno_training = unique(geno_training)
      
      geno_test =  merge(split[[2]], geno, by = 'geno_ID', all.x = T)[, col_to_keep]
      geno_test = unique(geno_test)
      
      
      rec <- recipe(geno_ID ~ . ,
                    data = geno_training) %>%
        update_role(geno_ID, new_role = 'outcome') %>%
        step_nzv(all_predictors()) %>%
        step_pca(
          all_predictors(),
          num_comp = num_pcs,
          options = list(center = T, scale. = T)
        )
      
      norm_obj <- prep(rec, training = geno_training)
      
      training_pca <- bake(norm_obj, geno_training)
      test_pca <- bake(norm_obj, geno_test)
      
      test_pca$geno_ID = geno_test$geno_ID
      
      pc_values <- rbind(training_pca, test_pca)
      pc_values <- unique(pc_values)
      pc_values <- as.data.frame(pc_values)
      
      
      return(pc_values)
      
    }
    
    pca_all_splits = lapply(splits,
                            apply_pca)
    
    
    
    splits = lapply(
      seq_along(splits),
      FUN = function(i) {
        append(splits[[i]], list(pca_all_splits[[i]]))
      }
    )
    
    
  }
  
  
  # Create the genotype matrix with SNP covariates selected if these shoud be
  # added.
  
  if (use_selected_markers == T) {
    SNPs = geno[, colnames(METData$geno) %in% METData$selected_markers]
    SNPs$geno_ID = rownames(SNPs)
  }
  
  
  ###############################
  ###############################
  
  
  ## PROCESSING AND SELECTING PREDICTORS FOR FITTING THE MODEL ##
  
  processing_all_splits = lapply(splits,
                                 function(x){processing_train_test_split(split = x,
                                                                         list_env_predictors = list_env_predictors,
                                                                         include_env_predictors = include_env_predictors,
                                                                         lat_lon_included = lat_lon_included,
                                                                         year_included = year_included,
                                                                         use_selected_markers = use_selected_markers)})
  
  ##  FITTING ALL TRAIN/TEST SPLITS OF THE EXTERNAL CV SCHEME
  
  fitting_all_splits = lapply(processing_all_splits,
                              function(x){fitting_train_test_split(split = x,
                                                                   prediction_method = method)})
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Select trait from the phenotypic file
  
  list_predictors_geno <- colnames(geno)
  
  pheno <-
    pheno[, c(trait, list_predictors_geno, list_env_predictors)]
  
  # Process the training and test sets
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Create the cross-validation random splits for hyperparameter optimization
  
  cv_splits <-
    rsample::vfold_cv(pheno, folds = nb_folds_cv, repeats = reps)
  
  # Define the type of model to tune
  
  mod <- parsnip::linear_reg(penalty = tune(),
                             mixture = tune()) %>%
    parsnip::set_engine("glmnet")
  
  
  # Define the predictors and outcome variables.
  # Remove predictors with null varaiance.
  # Center and scale the variables (standardization)
  
  rec <- recipe( ~ ., data = pheno) %>%
    update_role(all_of(trait), new_role = 'outcome') %>%
    update_role(all_of(list_predictors), new_role = 'predictor') %>%
    step_nzv(all_numeric()) %>%
    step_normalize(all_numeric())
  
  
}


#@reps:number of k - folds random partitions
