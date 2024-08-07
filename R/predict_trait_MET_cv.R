#' Cross-validation procedure for phenotypic prediction of crop varieties.
#'
#' @description
#' Implement trait prediction based on SNP and environmental data
#' with selection of prediction methods among Machine Learning approaches.
#' This function should be used to assess the predictive ability according to
#' a cross-validation scheme determined by the user.
#'
#' @param METData \code{list} An object created by the initial function of the
#' package [create_METData()].
#' @param trait \code{character} Name of the trait to predict.
#' @param prediction_method \code{character} specifying the predictive model to 
#' use.
#' Options are currently `xgb_reg_1` (gradient boosted trees), `xgb_reg_2` ,
#' `xgb_reg_3`, `DL_reg_1` (multilayer perceptrons), `DL_reg_2`, `DL_reg_3`,
#' `stacking_reg_1` (stacked models), `stacking_reg_2`, `stacking_reg_3`,
#' `rf_reg_1`, `rf_reg_2`, `rf_reg_3`.
#'
#' @param use_selected_markers A \code{Logical} indicating whether to use a
#' subset of markers  identified via single-environment GWAS or based on the
#' table of marker effects obtained via Elastic Net as predictor variables,
#' when main genetic effects are modeled with principal components. \cr
#' If `use_selected_markers` is `TRUE`, and if `list_selected_markers_manual`
#' is `NULL`, then the [select_markers()] function will be called in the
#' pipeline.
#' \strong{For more details, see [select_markers()]}
#'
#' @param year_included \code{logical} indicates if year factor should be used
#' as predictor variable. Default is `FALSE`.
#' 
#' @param location_included \code{logical} indicates if location factor should 
#'  be used as predictor variable. Default is `FALSE`.
#' 
#' @param lat_lon_included \code{logical} indicates if longitude and latitude 
#' data should be used as numeric predictors. Default is `FALSE`.
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
#'   of the environmental predictors which should be used in predictions.
#'   \strong{By default `NULL`: all environmental predictors included in the
#'   env_data table of the `METData` object will be used.}
#'
#' @param save_splits A \code{Logical} to indicate if the train/test splits
#'   should be saved.
#'
#' @param seed \code{integer} Seed value. Default is `NULL`. By default, a
#'   random seed will be generated.
#'
#' @param save_processing a \code{logical} indicating whether the processing
#'   steps obtained from the [get_splits_processed_with_method()] functions 
#'   should be saved in a .RDS object. Default is `FALSE`.
#'
#' @param path_folder a \code{character} indicating the full path where the .RDS
#'   object and plots generated during the analysis should be saved (do not use
#'   a Slash after the name of the last folder). Default is `NULL`.
#'
#' @param save_model a \code{logical} indicating Logical indicating whether the
#'   fitted model for each training-test partition should be saved. Default is
#'   FALSE. Note that some models (e.g. stacked models) can require a large
#'   memory.
#'
#' @param ... Arguments passed to the [get_splits_processed_with_method()]
#' function.
#'
#' @return A `list` object of class `met_cv` with the following items:
#'   \describe{
#'     \item{list_results_cv}{\code{list} of `res_fitted_split` elements.
#'     The length of this list corresponds to the number of training/test set
#'     partitions.}
#'     \item{seed_used}{\code{integer} Seed used to generate the
#'     cross-validation splits.}
#'     \item{cv_type}{\code{integer} Seed used to generate the
#'     cross-validation splits.}
#'     }
#'
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @export
#'
#'
predict_trait_MET_cv <- function(METData,
                                 trait,
                                 prediction_method,
                                 lat_lon_included = F,
                                 year_included = F,
                                 location_included = T,
                                 cv_type = "cv0",
                                 cv0_type = "leave-one-environment-out",
                                 nb_folds_cv1 = 5,
                                 repeats_cv1 = 50,
                                 nb_folds_cv2 = 5,
                                 repeats_cv2 = 50,
                                 nb_folds_cv00 = 5,
                                 include_env_predictors = T,
                                 list_env_predictors = NULL,
                                 use_selected_markers = F,
                                 list_selected_markers_manual = NULL,
                                 seed = NULL,
                                 save_splits = F,
                                 save_processing = F,
                                 path_folder,
                                 save_model = F,
                                 ...) {
  # Check the path_folder: create if does not exist
  path_folder <-
    file.path(path_folder,
              paste0(trait,
                     "_",
                     prediction_method,
                     "_",
                     cv_type))
  
  if (!dir.exists(path_folder)) {
    dir.create(file.path(path_folder), recursive = T)
  }
  
  
  # Test trait given by user
  if (is.null(trait)) {
    stop("Please give the name of the trait")
  }
  if (all(is.na(METData$pheno[, trait]))) {
    stop("Only NA values for this trait. Please check data or select",
         "another trait.")
  }
  
  # Define geno data
  
  geno <- as.data.frame(METData$geno)
  
  
  
  # Genotype matrix with SNP covariates selected if these should be added
  # as specific additional covariates (in addition to the main genetic effects
  # captured by the principal components).
  if (!use_selected_markers &
      length(list_selected_markers_manual) == 0 &
      prediction_method %in% c("stacking_reg_2")) {
    stop(
      cat(
        "stacking_reg_2 prediction method requires a subset of ",
        "marker variables. They can be provided via the argument ",
        "list_selected_markers_manual, or determined via the package",
        "using use_selected_markers = TRUE with Elastic Net or ",
        "GWAS approach."
      )
    )
  }
  
  if (use_selected_markers &
      length(list_selected_markers_manual) > 0) {
    SNPs <- as.data.frame(geno[, colnames(geno) %in% 
                                 list_selected_markers_manual])
    SNPs$geno_ID <- row.names(SNPs)
  } else if (use_selected_markers &
             length(list_selected_markers_manual) == 0) {
    list_selected_markers <- select_markers(
      METData = METData,
      trait = trait,
      path_save_res = file.path(path_folder, "GWAS"),
      ...
    )
    SNPs <- as.data.frame(geno[, colnames(geno) %in% list_selected_markers])
    SNPs$geno_ID <- row.names(SNPs)
    
  } else{
    cat("No specific SNP covariates will be used in analyses.\n")
  }
  
  
  # Check METData$ECs_computed to see if ECs were correctly downloaded via the
  # package when these are required by the user.
  
  if (include_env_predictors &
      is.null(METData$env_data)) {
    stop(
      cat(
        "No environmental covariates found in METData$env_data. Please set the",
        "argument compute_climatic_ECs to TRUE, or provide an environmental ",
        "data.frame"
      )
    )
  }
  # If no specific list of environmental predictors provided, all of the
  # environmental predictors present in METData$env_data are used as predictors.
  
  if (include_env_predictors &
      !is.null(METData$env_data)) {
    if (is.null(list_env_predictors)){
      list_env_predictors <- colnames(METData$env_data)[colnames(METData$env_data) %notin%
                                                          c("IDenv", "year", "location", "longitude", "latitude")]
    }
    
  }
  
  env_predictors <- METData$env_data
  info_environments <- METData$info_environments
  
  # Define the type of location information to use: categorical or lon-lat 
  # predictors
  
  if (location_included & !lat_lon_included){
    type_location_info <- "location_factor"
  }
  if (!location_included & lat_lon_included){
    type_location_info <- "lon_lat_numeric"
  }
  if (location_included & lat_lon_included){
    stop(
      cat(
        "Choose either the name of the location or its numeric coordinates",
        "to be included as predictors in the model.\n"
      )
    )
    
    
  }
  if (!location_included & !lat_lon_included){
    type_location_info <- "no_location_information"
    
  }
  
  # Select phenotypic data for the trait under study and remove NA in phenotypes
  
  pheno <-
    METData$pheno[, c("geno_ID", 
                      "year", 
                      "location", 
                      "IDenv", 
                      trait)][complete.cases(METData$pheno[, c("geno_ID", 
                                                               "year", 
                                                               "location", 
                                                               "IDenv", 
                                                               trait)]), ]
  
  # Create cross-validation random splits according to the type of selected CV
  
  # Generate a seed
  if (is.null(seed)) {
    seed_generated <- sample(size = 1, 1:2 ^ 15)
  } else{
    seed_generated <- seed
  }
  
  if (cv_type == "cv1") {
    cat(nb_folds_cv1,
        "-folds CV will be used for ",
        repeats_cv1,
        " repeats\n",
        sep = "")
    splits <-
      predict_cv1(
        pheno_data = pheno,
        nb_folds = nb_folds_cv1,
        reps = repeats_cv1,
        seed = seed_generated
      )
  }
  
  if (cv_type == "cv2") {
    cat(nb_folds_cv2,
        "-folds CV will be used for ",
        repeats_cv2,
        " repeats\n",
        sep = "")
    splits <-
      predict_cv2(
        pheno_data = pheno,
        nb_folds = nb_folds_cv2,
        reps = repeats_cv2,
        seed = seed_generated
      )
  }
  
  
  
  if (cv_type == "cv0") {
    splits <-
      predict_cv0(pheno_data = pheno,
                  cv0_type = cv0_type)
  }
  
  if (cv_type == "cv00") {
    splits <-
      predict_cv00(pheno_data = pheno,
                   cv0_type = cv0_type,
                   nb_folds = nb_folds_cv00)
  }
  
  if (cv_type == "cv00_5folds") {
    
    splits <-
      predict_cv00_5folds(pheno_data = pheno,
                          cv0_type = cv0_type,
                          nb_folds_cv00 = nb_folds_cv00)
  }
  
  if (save_splits) {
    saveRDS(splits,
            file = file.path(path_folder, "/splits.RDS"))
  }
  ###############################
  ###############################
  
  ## PROCESSING AND SELECTING PREDICTORS FOR FITTING THE MODEL ##
  
  
  checkmate::assert_class(splits,
                          "cv_object")
  
  checkmate::assert_choice(
    prediction_method,
    choices = c(
      "xgb_reg_1",
      "xgb_reg_2",
      "xgb_reg_3",
      "rf_reg_1",
      "rf_reg_2",
      "rf_reg_3",
      "DL_reg_1",
      "DL_reg_2",
      "DL_reg_3",
      "stacking_reg_1",
      "stacking_reg_2",
      "stacking_reg_3"
    )
  )
  
  processing_all_splits <-
    get_splits_processed_with_method(
      splits = splits,
      prediction_method = prediction_method,
      trait = trait,
      geno = geno,
      env_predictors = env_predictors,
      info_environments = info_environments,
      use_selected_markers = use_selected_markers,
      SNPs = SNPs,
      list_env_predictors = list_env_predictors,
      include_env_predictors = include_env_predictors,
      type_location_info = type_location_info,
      year_included = year_included,
      ...
    )
  
  
  if (save_processing) {
    saveRDS(processing_all_splits,
            file = file.path(path_folder, "/recipes_processing_splits.RDS"))
  }
  
  ###############################
  ###############################
  
  ##  FITTING ALL TRAINING SETS AND PREDICTING EACH TEST FOR EACH 
  ## SPLIT ELEMENT  ##
  
  fitting_all_splits <- list()
  length(fitting_all_splits) <- length(processing_all_splits)
  optional_args <- list(...)
  optional_args$seed <- seed_generated
  optional_args$path_folder <- path_folder
  optional_args$save_model <- save_model
  
  for (i in seq_len(length(fitting_all_splits))) {
    optional_args$object <- processing_all_splits[[i]]
    fitting_all_splits[[i]] <-
      do.call(fit_cv_split, args = optional_args)
  }
  #############################################
  ## SAVE RESULTS ALONG WITH THE SEED USED ####
  if (cv_type == "cv0") {
    cv_type <- paste0(cv_type, "_", cv0_type)
  } else{
    cv_type <- cv_type
  }
  
  met_cv <-
    list(
      "list_results_cv" = fitting_all_splits,
      "seed_used" = seed_generated,
      "cv_type" = cv_type
    )
  
  class(met_cv) <- c("list", "met_cv")
  
  saveRDS(met_cv,
          file = file.path(path_folder, "/met_cv.RDS"))
  
  ###############################
  ###############################
  
  ## VISUALIZATION OF THE PREDICTIVE ABILTIES AND OF 
  ## VARIABLE IMPORTANCE ACCORDING TO THE SELECTED CV SCHEME 
  
  plot_res <- plot_results_cv(
    fitting_all_splits = fitting_all_splits,
    trait = trait,
    info_environments = info_environments,
    cv_type = cv_type,
    cv0_type = cv0_type,
    path_folder = path_folder,
    nb_folds_cv1 = nb_folds_cv1,
    repeats_cv1 = repeats_cv1,
    nb_folds_cv2 = nb_folds_cv2,
    repeats_cv2 = nb_folds_cv2
  )
  
  return(met_cv)
}