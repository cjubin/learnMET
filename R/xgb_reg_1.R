#' Processing of a split object to get data ready to be used and fitted with
#' a `xgb_reg_1` (gradient boosted tree) regression model.
#'
#' @description
#' The function processes a split object (training + test sets), according to
#' the configuration set by the user. For instance, genomic information is
#' incorporated according to the option set by the user.  A list of specific
#' environmental covariables to use can be provided.\cr
#'
#' A recipe is created using the package `recipes`, to specify additional
#' preprocessing steps, such as standardization based on the training set, with
#' same transformations used on the test set. Variables with null variance are
#' removed. If year effect is included, it is converted to dummy variables.\cr
#' Further fitting on the training set with a gradient boosting model (see
#' function [fit_cv_split.xgb_reg_1()])).
#'
#'
#' @param split an object of class `split`. 
#'   A `split` object contains a training and test elements.
#'
#' @param trait \code{character} Name of the trait to predict. An ordinal trait
#'   should be encoded as `integer`.
#'
#' @param geno \code{data.frame} It corresponds to a `geno` element
#'   within an object of class `METData`.
#'
#' @param env_predictors \code{data.frame} It corresponds to the `env_data`
#'   element within an object of class `METData`.
#'
#' @param info_environments \code{data.frame} It corresponds to the
#'   `info_environments` element within an object of class `METData`.
#'
#' @param use_selected_markers A \code{Logical} indicating whether to use a
#'   subset of markers  identified via single-environment GWAS or based on the
#'   table of marker effects obtained via Elastic Net as predictor variables,
#'   when main genetic effects are modeled with principal components. \cr
#'   If `use_selected_markers` is `TRUE`, the `SNPs` argument should be
#'    provided.
#'    \strong{For more details, see [select_markers()]}
#'
#' @param SNPs A \code{data.frame} with the genotype matrix (individuals in rows
#'   and selected markers in columns) for SNPs selected via the
#'   [select_markers()] function.
#'   \strong{Optional argument, can remain as `NULL` if no single markers should
#'   be incorporated as predictor variables in analyses based on PCA
#'   decomposition.}
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
#' @param type_location_info \code{logical} indicates how the location 
#'   information is encoded, among `location_factor`, `lon_lat_numeric`, 
#'   `no_location_information`. Default is `location_factor`, meaning that the 
#'   variable location is used as a categorical variable in the model (encoded 
#'   as dummy variable).
#'
#' @param year_included \code{logical} indicates if year factor should be used
#'   as predictor variable. Default is `FALSE`.
#'
#' @return A `list` object of class `xgb_reg_1` with the following items:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set after partial processing}
#'     \item{test}{\code{data.frame} Test set after partial processing}
#'     \item{rec}{A \code{recipe} object, specifying the remaining processing
#'     steps which are implemented when a model is fitted on the training set
#'     with a recipe.}
#' }
#'
#' @references
#' \insertRef{wickham2019welcome}{learnMET}
#' \insertRef{tidymodels}{learnMET}
#'
#' @name xgb_reg_1
#' @export
new_xgb_reg_1 <- function(split = NULL,
                          trait = NULL,
                          geno = NULL,
                          env_predictors = NULL,
                          info_environments = NULL,
                          use_selected_markers = F,
                          SNPs = NULL,
                          include_env_predictors = T,
                          list_env_predictors = NULL,
                          type_location_info = "location_factor",
                          year_included = F,
                          ...) {
  if (class(split) != 'split') {
    stop('Class of x should be "split".')
  }
  if (class(split[['training']][, trait]) %in% c('integer')) {
    split[['training']][, trait] <-
      as.numeric(split[['training']][, trait])
    split[['test']][, trait] <- as.numeric(split[['test']][, trait])
    
  }
  
  
  ## GENOTYPIC DATA ##
  
  # Use of genotypic data: use of PCs derived from additive geno matrix #
  cat('Processing: PCA transformation on the scaled marker dataset\n')
  
  
  
  pca_geno = apply_pca(split = split,
                       geno = geno)
  
  training = pca_geno[[1]]
  test = pca_geno[[2]]
  cat('Processing: PCA transformation done\n')
  
  
  
  # Add specific SNP covariates if they should be used
  
  if (use_selected_markers) {
    training <- merge(training, SNPs, by = 'geno_ID', all.x = T)
    test <- merge(test, SNPs, by = 'geno_ID', all.x = T)
    
  }
  
  
  ## ENVIRONMENTAL DATA ##
  # Add the environmental data
  
  if (include_env_predictors &
      !is.null(list_env_predictors)) {
    training <-
      merge(training,
            env_predictors[, c('IDenv', list_env_predictors)],
            by = 'IDenv',
            all.x = T)
    test <-
      merge(test,
            env_predictors[, c('IDenv', list_env_predictors)],
            by = 'IDenv',
            all.x = T)
    
  }
  
  
  if (type_location_info == "lon_lat_numeric" &
      year_included &
      length(unique(as.character(training$year))) > 1) {
    # Add longitude/latitude data for each train & test split
    
    training <-
      merge(training,
            info_environments[, c('IDenv', 'longitude', 'latitude')],
            by = 'IDenv',
            all.x = T)
    test <-
      merge(test,
            info_environments[, c('IDenv', 'longitude', 'latitude')],
            by = 'IDenv',
            all.x = T)
    
    
    
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_novel(year,location,geno_ID,IDenv) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(geno_ID) %>%
      recipes::update_role(-trait, -IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, keep_original_cols = F, one_hot = TRUE) %>%
      recipes::step_nzv(recipes::all_predictors(), -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes(),-tidyselect::starts_with('PC'))
    
    
    
    
    
  }
  else if (type_location_info == "location_factor" &
           year_included &
           length(unique(as.character(training$year))) > 1) {
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_novel(year,location,geno_ID,IDenv) %>%
      recipes::step_rm(geno_ID) %>%
      recipes::update_role(-trait, -IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, keep_original_cols = F, one_hot = TRUE) %>%
      recipes::step_dummy(location, keep_original_cols = F, one_hot = TRUE) %>%
      recipes::step_zv(recipes::all_predictors(), -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes(),-tidyselect::starts_with('PC'))
    
    
    
  }
  
  else if (type_location_info == "location_factor" &
           !year_included) {
    
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_rm(geno_ID) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-trait, -IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(location, keep_original_cols = F, one_hot = TRUE) %>%
      recipes::step_zv(recipes::all_predictors()) %>%
      recipes::step_nzv(recipes::all_predictors(),
                        -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(
        recipes::all_numeric(),
        -recipes::all_outcomes(),
        -tidyselect::starts_with('PC')
      )
    
    
  }
  else if ((type_location_info == "lon_lat_numeric" &
            !year_included) | (type_location_info == "lon_lat_numeric" &
                               length(unique(as.character(training$year))) <
                               2)) {
    # Add longitude/latitude data for each train & test split
    
    training <-
      dplyr::left_join(training,
                       info_environments[, c('IDenv', 'longitude', 'latitude')],
                       by = 'IDenv')
    test <-
      dplyr::left_join(test,
                       info_environments[, c('IDenv', 'longitude', 'latitude')],
                       by = 'IDenv')
    
    
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(geno_ID) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-trait, -IDenv, new_role = 'predictor') %>%
      recipes::step_zv(recipes::all_predictors()) %>%
      recipes::step_nzv(recipes::all_predictors(),
                        -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(
        recipes::all_numeric(),
        -recipes::all_outcomes(),
        -tidyselect::starts_with('PC')
      )
    
    
    
    
  } else if (type_location_info == "no_location_information" &
             year_included &
             length(unique(as.character(training$year))) > 1) {
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_rm(geno_ID) %>%
      recipes::step_rm(location) %>%
      recipes::update_role(-trait, -IDenv,new_role = 'predictor') %>%
      recipes::step_dummy(year, keep_original_cols = F, one_hot = TRUE) %>%
      recipes::step_zv(recipes::all_predictors()) %>%
      recipes::step_nzv(recipes::all_predictors(),
                        -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(
        recipes::all_numeric(),
        -recipes::all_outcomes(),
        -tidyselect::starts_with('PC')
      )
    
    
  } 
  else{
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipes::recipe( ~ . ,
                            data = training) %>%
      recipes::update_role(trait, new_role = 'outcome') %>%
      recipes::update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(geno_ID) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-trait, -IDenv, new_role = 'predictor') %>%
      recipes::step_zv(recipes::all_predictors()) %>%
      recipes::step_nzv(recipes::all_predictors(),
                        -tidyselect::starts_with('PC')) %>%
      recipes::step_normalize(
        recipes::all_numeric(),
        -recipes::all_outcomes(),
        -tidyselect::starts_with('PC')
      )
    
    
  }
  cat(
    'Incorporating selected predictors & Data processing for one train/test split: Done!\n'
  )
  
  split_processed <- structure(list(
    "training" = training,
    "test" = test,
    "rec" = rec
  ), class = 'xgb_reg_1')
  
  
  
  
  return(split_processed)
  
  
}






#' @rdname xgb_reg_1
#' @aliases new_xgb_reg_1
#' @export
xgb_reg_1 <- function(split,
                      trait,
                      geno,
                      env_predictors,
                      info_environments,
                      use_selected_markers,
                      SNPs,
                      list_env_predictors,
                      include_env_predictors,
                      type_location_info,
                      year_included,
                      ...) {
  validate_xgb_reg_1(
    new_xgb_reg_1(
      split = split,
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
  )
}


#' @rdname xgb_reg_1
#' @aliases new_xgb_reg_1
#' @export

validate_xgb_reg_1 <- function(x, ...) {
  trait <-
    as.character(x[['rec']]$term_info[which(x[['rec']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'xgb_reg_1')
  
  checkmate::assert_names(names(x), must.include = c('training', 'test', 'rec'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  
  
  return(x)
}