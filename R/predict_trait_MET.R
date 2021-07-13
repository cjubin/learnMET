#' Cross-validation procedure for phenotypic prediction of crop varieties.
#'
#' @description
#' Implement trait prediction based on SNP and environmental data
#' with selection of prediction methods among Machine Learning approaches.
#'
#' This function should be used to assess the predictive ability according to
#' a cross-validation scheme determined by the user.
#'
#' @param METData \code{list} An object created by the initial function of the
#'   package create_METData().
#'
#' @param trait \code{character} Name of the trait to predict. An ordinal trait
#'   should be encoded as `integer`.
#'
#' @param method_processing \code{character} specifying the predictive model to use.
#'   Options are `xgb_reg` (gradient boosted trees),
#'   (stacking of support vector machines with LASSO as meta-learner).
#'
#' @param use_selected_markers A \code{Logical} indicating whether to use a
#'   subset of markers obtained from a previous step
#'   (see [function select_markers()]).
#'
#' @param geno_information indicating how the complete genotypic matrix should
#'   be used in predictions. Options are `SNPs` (all
#'   of the markers will be individually used), `PCs` (PCA will be applied on
#'   each genotype matrix for the training set for dimensionality reduction)
#'   or `PCs_G` (decomposition of the genomic relationship matrix via PCA -- not
#'   yet implemented).
#'
#' @param num_pcs \code{}. Default is 200.
#'
#' @param lat_lon_included \code{logical} indicates if longitude and latitude
#'   data should be used as numeric predictors. Default is `TRUE`.
#'
#' @param year_included \code{logical} indicates if year factor should be used
#'   as predictor variable. Default is `FALSE`.
#'
#' @param cv_type A \code{character} with one out of `cv0` (prediction of new
#'   environments), `cv00` (prediction of new genotypes in new environments),
#'   `cv1` (prediction of new genotypes) or `cv2` (prediction of incomplete
#'   field trials). Default is `cv0`.
#'
#' @param cv0_type A \code{character} with one out of
#'   `leave-one-environment-out`, `leave-one-site-out`,`leave-one-year-out`,
#'   `forward-prediction`. Default is `leave-one-environment-out`.
#'
#' @param nb_folds_cv1 A \code{numeric} Number of folds used in the CV1 scheme.
#'   Default is 5.
#'
#' @param repeats_cv1 A \code{numeric} Number of repeats in the CV1 scheme.
#'   Default is 50.
#'
#' @param nb_folds_cv2 A \code{numeric} Number of folds used in the CV2 scheme.
#'   Default is 5.
#'
#' @param repeats_cv2 A \code{numeric} Number of repeats in the CV2 scheme.
#'   Default is 50.
#'
#' @param include_env_predictors A \code{logical} indicating whether
#'   environmental covariates characterizing each environment should be used in
#'   predictions.
#'
#' @param list_env_predictors A \code{character} vector containing the names
#'   of the environmental predictors which should be used in predictions. By
#'   default `NULL`: all environmental predictors included in the env_data table
#'   of the `METData` object will be used.
#'
#' @param seed \code{integer} Seed value. Default is `NULL`. By default, a
#'   random seed will be generated.
#'
#' @param save_processing a \code{logical} indicating whether the processing
#'   steps obtained from the [processing_train_test_split()] or
#'   [processing_train_test_split_kernel()] functions should be saved in a .RDS
#'   object. Default is `FALSE`.
#'
#' @param path_folder a \code{character} indicating the full path where the .RDS
#'   object and plots generated during the analysis should be saved (do not use
#'   a Slash after the name of the last folder). Default is `NULL`.
#'
#' @param ... Arguments passed to the [processing_train_test_split()],
#'  [processing_train_test_split_kernel()], [reg_fitting_train_test_split()],
#'  [reg_fitting_train_test_split_kernel()] functions.
#'
#' @return A `list` object of class `met_cv` with the following items:
#'   \describe{
#'     \item{list_results_cv}{\code{list} of `res_fitted_split` elements.
#'     Detailed prediction results for each split of the
#'     data within each element of this list.}
#'     \item{seed_used}{\code{integer} Seed used to generate the
#'     cross-validation splits.}
#'     }
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'
predict_trait_MET <- function(METData_training,
                              METData_new,
                              trait,
                              method_processing = c('xgb_reg'),
                              use_selected_markers = F,
                              geno_information = c('PCs'),
                              num_pcs = 200,
                              lat_lon_included = F,
                              year_included = F,
                              include_env_predictors = T,
                              list_env_predictors = NULL,
                              seed = NULL,
                              save_processing = T,
                              path_folder,
                              vip_plot = TRUE,
                              ...) {
  
  
  
  # Check classes of METData
  checkmate::assert_class(METData_training, 'METData')
  checkmate::assert_class(METData_new, 'METData')
  
  # Check the path_folder: create if does not exist
  
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = T)
  }
  
  
  # Test trait given by user
  if (is.null(trait)) {
    stop('Please give the name of the trait')
  }
  
  # Define geno data 
  
  geno = METData_new$geno

  
  # Genotype matrix with SNP covariates selected if these should be added
  # as specific additional covariates (in addition to the main genetic effects).
  
  if (use_selected_markers == T &
      length(METData$selected_markers) > 0) {
    SNPs = as.data.frame(geno[, colnames(geno) %in% METData$selected_markers])
    SNPs$geno_ID = row.names(SNPs)
  } else if (use_selected_markers == T &
             length(METData$selected_markers) == 0) {
    cat(
      paste(
        'SNP covariates required to be used but no list of selected markers',
        'available.\nPlease use select_markers() and use the output of this',
        'function to run predict_trait_MET_cv() with selected SNP covariates.'
      )
    )
  } else{
    cat('No specific additional SNP covariates will be used in analyses.')
  }
  
  
  # Check METData_training$env_data and METData_new$env_data
  
  if (include_env_predictors &
      is.null(METData_training$env_data) & !METData_training$compute_ECs) {
    stop(
      'No environmental covariates found in METData_training$env_data. Please', 
      'set the argument "compute_ECs" to TRUE or provide environmental data.'
    )
  }
  
  if (include_env_predictors &
      is.null(METData_new$env_data) & !METData_new$compute_ECs) {
    stop(
      'No environmental covariates found in METData_training$env_data. Please', 
      'set the argument "compute_ECs" to TRUE or provide environmental data.'
    )
  }
  
  # Use only common columns in the env_data objects from training data
  # and from data to predict --> set of environmental predictors.
  common_cols <-
    intersect(colnames(METData_training$env_data), colnames(METData_new$env_data))
  METData_training$env_data <- METData_training$env_data[, common_cols]
  METData_new$env_data <- METData_new$env_data[, common_cols]
  
  
  # If no specific list of environmental predictors provided, all of the
  # environmental predictors present in METData_training$env_data are used as predictors.
  if (is.null(list_env_predictors) &
      include_env_predictors &
      nrow(METData_training$env_data) > 0 &  nrow(METData_new$env_data) > 0) {
    list_env_predictors = colnames(METData_training$env_data)[colnames(METData_training$env_data) %notin%
                                                                c('IDenv', 'year', 'location', 'longitude', 'latitude')]
    
    
  }
  
  env_predictors = rbind(METData_new$env_data,METData_training$env_data)
  
  
  # Select phenotypic data for the trait under study and remove NA in phenotypes
  
  pheno_training = METData_training$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)][complete.cases(METData_training$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)]),]
  pheno_new[,trait]=NA

  # Create cross-validation random splits according to the type of selected CV
  
  # Generate a seed
  if (is.null(seed)) {
    seed_generated <- sample(size = 1, 1:2 ^ 15)
  }
  else{
    seed_generated <- seed
  }
  
  split <- list(METData_training$pheno, METData_new$pheno)
  names(split) <- c('training', 'test')
  class(split) <- 'split'
  
  
  ###############################
  ###############################
  
  ## PROCESSING AND SELECTING PREDICTORS FOR FITTING THE MODEL ##
  #names_selected_SNPs <- colnames(SNPs)[colnames(SNPs) %notin% 'geno_ID']
  
  processing_tr_te_sets <-
    get_splits_processed_with_method(
      splits = split,
      method_processing = method_processing,
      trait = trait,
      geno_data = geno,
      env_predictors = env_predictors,
      info_environments = METData$info_environments,
      unique_EC_by_geno = METData$unique_EC_by_geno,
      geno_information = geno_information,
      use_selected_markers = use_selected_markers,
      SNPs = SNPs,
      list_env_predictors = list_env_predictors,
      include_env_predictors = include_env_predictors,
      lat_lon_included = lat_lon_included,
      year_included = year_included
    )
  
  
  
  if (save_processing) {
    saveRDS(
      processing_tr_te_sets,
      file = file.path(
        path_folder,
        '/recipes_processing_completeTrSetNewData.RDS'
      )
    )
  }
  
  ###############################
  ###############################
  
  ##  Fitting the complete METData training object and predicting the test set ##
  
  
  fit_and_predictions = fit_predict_model(object = processing_tr_te_sets,
                             seed = seed_generated,
                             path_folder = path_folder,
                             ...)
  
  
  ###############################
  ###############################
  
  
  ## VISUALIZATION OF THE VARIABLE IMPORTANCE FROM THE TRAINING SET ##
  
  if (vip_plot) {
    plot_vip <- plot_results_vip(
      fitting_all_splits = fit_and_predictions,
      method_processing = method_processing,
      path_folder = path_folder
    )
    
  }
  
  
  ## RETURNING RESULTS ALONG WITH THE SEED USED
  
  met_pred <-
    list(
      'list_results' = fit_and_predictions,
      'seed_used' = seed_generated
    )
  
  class(met_pred) <- c('list', 'met_pred')
  
  return(met_pred)
  
  
}
