#' Processing of a split object to get data ready to be used and fitted with
#' a `stacking_reg_1` (stacking of SVM models) regression model.
#'
#' @description
#' The function processes a split object (training + test sets), according to
#' the configuration set by the user. For instance, genomic information is 
#' incorporated according to the option set by the user.  A list of specific
#' environmental covariables to use can be provided.\cr
#'
#' Multiple recipes are created using the package `recipes` according to the
#' data source (genomic, environmental)
#' These recipes specify additional preprocessing steps, such as standardization
#' based on the training set, with same transformations used on the test set.
#' Variables with null variance are removed. If year effect is included, it is
#' converted to dummy variables. \cr Each recipe (with G and E features) will be
#' subsequently fitted with a support vector regression model and predictions
#' from each model will be combined (see function 
#' [fit_cv_split.stacking_reg_1()]).
#' \cr
#'
#' @param an object of class `split`, which is a subelement of the output of the
#'   [predict_cv00()], [predict_cv0()], [predict_cv1()] and [predict_cv2()]
#'   functions. A `split` object contains a training and test elements.
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
#'   subset of markers  identified via single-environment GWAS or based on the
#'   table of marker effects obtained via Elastic Net as predictor variables,
#'   when main genetic effects are modeled with principal components. \cr
#'   If `use_selected_markers` is `TRUE`, the `SNPs` argument should be
#'   provided.
#'   \strong{For more details, see [select_markers()]}
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
#' @param lat_lon_included \code{logical} indicates if longitude and latitude
#'   data should be used as numeric predictors. Default is `FALSE`.
#'
#' @param year_included \code{logical} indicates if year factor should be used
#'   as predictor variable. Default is `FALSE`.
#'
#' @return A `list` object of class `stacking_reg_1` with the following items:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set after partial processing}
#'     \item{test}{\code{data.frame} Test set after partial processing}
#'     \item{rec_G}{A \code{recipe} object, specifying the remaining processing
#'     steps which are implemented when a model is fitted on the training set
#'     with a recipe. Data used are predictors corresponding to genomic data.}
#'     \item{rec_E}{A \code{recipe} object, specifying the remaining processing
#'     steps which are implemented when a model is fitted on the training set
#'     with a recipe. Data used are predictors corresponding to enviornmental
#'     predictors.}
#'   }
#'
#' @references
#' \insertRef{wickham2019welcome}{learnMET}
#' \insertRef{tidymodels}{learnMET}
#'
#' @name stacking_reg_1
#' @export
new_stacking_reg_1 <- function(split = NULL,
                               trait = NULL,
                               geno = NULL,
                               env_predictors = NULL,
                               info_environments = NULL,
                               use_selected_markers = F,
                               SNPs = NULL,
                               include_env_predictors = T,
                               list_env_predictors = NULL,
                               lat_lon_included = F,
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
  
  
  geno$geno_ID = row.names(geno)
  
  ## SNPs DATA ##
  # Add the genotype data
  
  # Merge in same data.frame pheno and geno data for each train & test split
  
  training <-
    merge(split[[1]], geno, by = 'geno_ID', all.x = T)
  
  test <-
    merge(split[[2]], geno, by = 'geno_ID', all.x = T)
  
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
  
  
  ## ENVIRONMENTAL-BASED KERNEL ##
  
  if (lat_lon_included &
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
    
    rec_E <- recipes::recipe( ~ . ,
                              data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno))) %>%
      recipes::step_rm(location) %>%
      recipes::update_role(-tidyselect::all_of(trait), -IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes())
    
    
    
    
    
  } else if (!lat_lon_included &
             year_included &
             length(unique(as.character(training$year))) > 1) {
    # Create recipe to define the processing of the training & test set.
    
    rec_E <- recipes::recipe( ~ . ,
                              data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno))) %>%
      recipes::step_rm(location) %>%
      recipes::update_role(-tidyselect::all_of(trait), -IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes())
    
    
    
    
  } else if ((lat_lon_included &
              !year_included) | (lat_lon_included &
                                 length(unique(as.character(training$year))) <
                                 2)) {
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
    
    rec_E <- recipes::recipe( ~ . ,
                              data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno))) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-tidyselect::all_of(trait), -IDenv, new_role = 'predictor') %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes())
    
    
    
    
  } else{
    # Create recipe to define the processing of the training & test set.
    
    rec_E <- recipes::recipe( ~ . ,
                              data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno))) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-tidyselect::all_of(trait), -IDenv, new_role = 'predictor') %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes())
    
    
    
    
  }
  
  cat('Processing: recipe for the environmental-based kernel created!\n')
  
  ## GENOMIC BASED KERNEL ##
  
  rec_G <- recipes::recipe( ~ . ,
                            data = training) %>%
    recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
    recipes::update_role(IDenv, new_role = "id variable") %>%
    recipes::update_role(geno_ID, new_role = "id variable") %>%
    recipes::step_rm(location) %>%
    recipes::step_rm(year) %>%
    recipes::step_rm(all_of(list_env_predictors)) %>%
    recipes::update_role(-tidyselect::all_of(trait), -IDenv, -geno_ID, new_role = 'predictor') %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(),
    #                   skip = TRUE,
    #                   threshold = 0.95) %>%
    recipes::step_normalize(recipes::all_numeric(), -recipes::all_outcomes())
  
  
  
  cat('Processing: recipe for the genomic-based kernel created!\n')
  
  
  
  
  split_processed <- structure(
    list(
      'training' = training,
      'test' = test,
      'rec_G' = rec_G,
      'rec_E' = rec_E
    ),
    class = 'stacking_reg_1'
  )
  
  
  
  return(split_processed)
  
  
}











#' @rdname stacking_reg_1
#' @aliases new_stacking_reg_1
#' @export
stacking_reg_1 <- function(split,
                           trait,
                           geno,
                           env_predictors,
                           info_environments,
                           use_selected_markers,
                           SNPs,
                           list_env_predictors,
                           include_env_predictors,
                           lat_lon_included,
                           year_included,
                           ...) {
  validate_stacking_reg_1(
    new_stacking_reg_1(
      split = split,
      trait = trait,
      geno = geno,
      env_predictors = env_predictors,
      info_environments = info_environments,
      use_selected_markers = use_selected_markers,
      SNPs = SNPs,
      list_env_predictors = list_env_predictors,
      include_env_predictors = include_env_predictors,
      lat_lon_included = lat_lon_included,
      year_included = year_included,
      ...
    )
  )
}


#' @rdname stacking_reg_1
#' @aliases new_stacking_reg_1
#' @export
validate_stacking_reg_1 <- function(x, ...) {
  trait <-
    as.character(x[['rec_G']]$term_info[which(x[['rec_G']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'stacking_reg_1')
  
  checkmate::assert_names(names(x),
                          must.include = c('training', 'test', 'rec_G', 'rec_E'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  
  return(x)
}
