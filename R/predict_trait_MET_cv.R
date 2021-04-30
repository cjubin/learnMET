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
  
  
  ## Use of genotypic data: dimensionality reduction specified via parameter geno_information ##
  
  if (geno_information == 'PCs') {
    # Processing of PCs: apply transformations calculated on the training set,
    # on test set.
    
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
  
  ##  FITTING ALL TRAIN/TEST SPLITS OF THE EXTERNAL CV SCHEME ##
  
  fitting_all_splits = lapply(processing_all_splits,
                              function(x){fitting_train_test_split(split = x,
                                                                   prediction_method = method)})
  
  
  ## VISUALIZATION OF THE PREDICTIVE ABILTIES ACCORDING TO THE SELECTED CV SCHEME ##
  
  if (cv_type=='cv0'){
    if(cv0_type == 'leave-one-environment-out'){
      
    }
  }
  
  
}



