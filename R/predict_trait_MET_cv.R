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
predict_trait_MET_cv <- function(METData,
                                 trait,
                                 method_processing = c('xgb_reg'),
                                 use_selected_markers = T,
                                 geno_information = c('PCs'),
                                 num_pcs = 200,
                                 lat_lon_included = F,
                                 year_included = F,
                                 cv_type = c('cv0'),
                                 cv0_type = 'leave-one-environment-out',
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 include_env_predictors = T,
                                 list_env_predictors = NULL,
                                 seed = NULL,
                                 save_processing = F,
                                 path_folder = NULL,
                                 vip_plot = TRUE,
                                 ...) {
  # Check the path_folder: create if does not exist
  if (!dir.exists(path_folder)) {
    dir.create(path_folder, recursive = T)
  }
  
  
  # Test trait given by user
  if (is.null(trait)) {
    stop('Please give the name of the trait')
  }
  
  geno = METData$geno
  
  
  
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
  
  
  # Check METData$ECs_computed to see if ECs were correctly downloaded via the
  # package when these are required by the user.
  
  if (include_env_predictors &
      METData$compute_ECs & "ECs_computed" %notin% names(METData)) {
    stop(
      paste(
        'The weather-based covariates were not computed. Please use the',
        'function get_ECs() to obtain environmental predictors.'
      )
    )
  }
  
  if (include_env_predictors &
      is.null(METData$env_data)) {
    stop(
      'No environmental covariates found in METData$env_data. Please set the
      argument "compute_ECs" to TRUE when using create_METData(), and then run
      function get_ECs() to obtain environmental predictors based on weather
      data retrieved from NASA-POWER.'
    )
  }
  # If no specific list of environmental predictors provided, all of the
  # environmental predictors present in METData$env_data are used as predictors.
  if (is.null(list_env_predictors) &
      include_env_predictors & nrow(METData$env_data) > 0) {
    list_env_predictors = colnames(METData$env_data)[colnames(METData$env_data) %notin%
                                                       c('IDenv', 'year', 'location', 'longitude', 'latitude')]
    
    
  }
  
  env_predictors = METData$env_data
  
  
  # Select phenotypic data for the trait under study and remove NA in phenotypes
  
  pheno = METData$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)][complete.cases(METData$pheno[, c("geno_ID", "year" , "location", "IDenv", trait)]), ]
  
  
  # Create cross-validation random splits according to the type of selected CV
  
  # Generate a seed
  if (is.null(seed)) {
    seed_generated <- sample(size = 1, 1:2 ^ 15)
  }
  else{
    seed_generated <- seed
  }
  
  if (cv_type == 'cv1') {
    splits <-
      predict_cv1(
        pheno_data = pheno,
        nb_folds = nb_folds_cv1,
        reps = repeats_cv1,
        seed = seed_generated
      )
  }
  
  if (cv_type == 'cv2') {
    splits <-
      predict_cv2(
        pheno_data = pheno,
        nb_folds = nb_folds_cv2,
        reps = repeats_cv2,
        seed = seed_generated
      )
  }
  
  
  
  if (cv_type == 'cv0') {
    splits <-
      predict_cv0(pheno_data = pheno,
                  type = cv0_type)
  }
  print(length(splits))
  ###############################
  ###############################
  
  ## PROCESSING AND SELECTING PREDICTORS FOR FITTING THE MODEL ##
  #names_selected_SNPs <- colnames(SNPs)[colnames(SNPs) %notin% 'geno_ID']
  
  processing_all_splits <-
    get_splits_processed_with_method(
      splits = splits,
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
    saveRDS(processing_all_splits,
            file = file.path(path_folder, '/recipes_processing_splits.RDS'))
  }
  
  ###############################
  ###############################
  
  ##  FITTING ALL TRAIN/TEST SPLITS  ##
  
  
  
  fitting_all_splits = lapply(processing_all_splits,
                              function(x, ...) {
                                fit_cv_split(object = x,
                                             seed = seed_generated,
                                             path_folder = path_folder,
                                             ...)
                              })
  
  ###############################
  ###############################
  
  
  ## VISUALIZATION OF THE PREDICTIVE ABILTIES ACCORDING TO THE SELECTED CV SCHEME ##
  
  plot_res <- plot_results_cv(
    fitting_all_splits = fitting_all_splits,
    info_environments = METData$info_environments,
    method_processing = method_processing,
    splits,
    cv_type = cv_type,
    cv0_type = cv0_type,
    path_folder = path_folder,
    nb_folds_cv1 = nb_folds_cv1,
    repeats_cv1 = repeats_cv1,
    nb_folds_cv2 = nb_folds_cv2,
    repeats_cv2 = nb_folds_cv2
  )
  
  
  ## VISUALIZATION OF THE VARIABLE IMPORTANCE ##
  if (vip_plot) {
    plot_vip <- plot_results_vip(
      fitting_all_splits = fitting_all_splits,
      info_environments = METData$info_environments,
      method_processing = method_processing,
      splits,
      cv_type = cv_type,
      cv0_type = cv0_type,
      path_folder = path_folder,
      nb_folds_cv1 = nb_folds_cv1,
      repeats_cv1 = repeats_cv1,
      nb_folds_cv2 = nb_folds_cv2,
      repeats_cv2 = nb_folds_cv2
    )
    
  }
  
  if (cv_type == 'cv0') {
    type_cv <- paste0(cv_type, '_', cv0_type)
  } else{
    type_cv <- cv_type
  }
  
  ## RETURNING RESULTS ALONG WITH THE SEED USED
  met_cv <-
    list(
      'list_results_cv' = fitting_all_splits,
      'seed_used' = seed_generated,
      'type_cv' = type_cv
    )
  
  class(met_cv) <- c('list', 'met_cv')
  
  return(met_cv)
  
  
}
