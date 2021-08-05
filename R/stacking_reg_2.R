#' Builds a recipe to process a split object (containing training
#' and test sets) according to the configuration set by the user and assign it
#' to a model that stacks two or three support vector regression models (SVM)
#' for subsequent model fitting using a S3 method dispatch.
#'
#' @description
#' The function processes genomic information according to the option set by the
#' user. Training and test datasets are subsetted on columns based on the
#' list of environmental variables to use.\cr
#'
#' Multiple recipes are created using the package `recipes` according to the
#' data source (genomic, environmental or first-order GxE interactions datasets)
#' These recipes specify additional preprocessing steps, such as standardization
#' based on the training set, with same transformations used on the test set.
#' Variables with null variance are removed. If year effect is included, it is
#' converted to dummy variables. \cr Each recipe (with G, E or GxE) will be
#' subsequently fitted with a support vector regression model and predictions
#' from each model will be combined (see function
#' [fit_cv_split.stacking_reg_2()]).
#' \cr
#'
#' @param an object of class `split`, which is a subelement of the output of the
#'   [predict_cv00()], [predict_cv0()], [predict_cv1()] and [predict_cv2()]
#'   functions. A `split` object contains a training and test elements.
#'
#' @param trait \code{character} Name of the trait to predict. An ordinal trait
#'   should be encoded as `integer`.
#'
#' @param geno_data \code{data.frame} It corresponds to a `geno` element
#'   within an object of class `METData`.
#'
#' @param env_predictors \code{data.frame} It corresponds to the `env_data`
#'   element within an object of class `METData`.
#'
#' @param info_environments \code{data.frame} It corresponds to the
#'   `info_environments` element within an object of class `METData`.
#'
#' @param geno_information A \code{character} indicating how the complete
#'   genotype matrix should be used in predictions. Options are `SNPs` (all
#'   of the markers will be individually used), `PCs` (PCA will be applied on
#'   each genotype matrix for the training set for dimensionality reduction)
#'   or `PCs_G` (decomposition of the genomic relationship matrix via eigen
#'   value decomposition).
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
#' @return A `list` object of class `stacking_reg_2` with the following items:
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
#'     \item{rec_GE}{A \code{recipe} object, specifying the remaining processing
#'     steps which are implemented when a model is fitted on the training set
#'     with a recipe. Data used are predictors corresponding to first-order
#'     GxE interactions.}
#'   }
#'
#' @references
#' \insertRef{wickham2019welcome}{learnMET}
#' \insertRef{tidymodels}{learnMET}
#'
#' @name stacking_reg_2
#' @export
new_stacking_reg_2 <- function(split = NULL,
                               trait = NULL,
                               geno_data = NULL,
                               env_predictors = NULL,
                               info_environments = NULL,
                               geno_information = 'SNPs',
                               use_selected_markers = T,
                               SNPs = NULL,
                               include_env_predictors = T,
                               list_env_predictors = NULL,
                               lat_lon_included = F,
                               year_included = F,
                               ...) {
  if (class(split) != 'split') {
    stop('Class of x should be "split".')
  }
  
  if (!use_selected_markers & is.null(SNPs)) {
    stop(paste(
      'Method "stacking_reg_2" uses a subset of SNPs along with ECs ',
      'in a kernel. To use this method, please set the argument ',
      'use_selected_markers to TRUE.'
    ))
  }
  
  if (class(split[['training']][, trait]) %in% c('integer')) {
    split[['training']][, trait] <-
      as.numeric(split[['training']][, trait])
    split[['test']][, trait] <- as.numeric(split[['test']][, trait])
    
  }
  
  
  geno_data$geno_ID = row.names(geno_data)
  
  ## SNPs DATA ##
  # Add the genotype data
  
  # Merge in same data.frame pheno and geno data for each train & test split
  
  training <-
    merge(split[[1]], geno_data, by = 'geno_ID', all.x = T)
  
  test <-
    merge(split[[2]], geno_data, by = 'geno_ID', all.x = T)
  
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
    
    rec_E <- recipes::recipe(~ . ,
                             data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno_data))) %>%
      recipes::step_rm(location) %>%
      recipes::update_role(-tidyselect::all_of(trait),-IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
    
    
    
    
    
  } else if (!lat_lon_included &
             year_included &
             length(unique(as.character(training$year))) > 1) {
    # Create recipe to define the processing of the training & test set.
    
    rec_E <- recipes::recipe(~ . ,
                             data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno_data))) %>%
      recipes::step_rm(location) %>%
      recipes::update_role(-tidyselect::all_of(trait),-IDenv, new_role = 'predictor') %>%
      recipes::step_dummy(year, preserve = F, one_hot = TRUE) %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
    
    
    
    
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
    
    rec_E <- recipes::recipe(~ . ,
                             data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno_data))) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-tidyselect::all_of(trait),-IDenv, new_role = 'predictor') %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
    
    
    
    
  } else{
    # Create recipe to define the processing of the training & test set.
    
    rec_E <- recipes::recipe(~ . ,
                             data = training) %>%
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_rm(all_of(colnames(geno_data))) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-tidyselect::all_of(trait),-IDenv, new_role = 'predictor') %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
    
    
    
    
  }
  
  cat('Processing: recipe for the environmental-based kernel created!\n')
  
  ## GENOMIC BASED KERNEL ##
  
  rec_G <- recipes::recipe(~ . ,
                           data = training) %>%
    recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
    recipes::update_role(IDenv, new_role = "id variable") %>%
    recipes::update_role(geno_ID, new_role = "id variable") %>%
    recipes::step_rm(location) %>%
    recipes::step_rm(year) %>%
    recipes::step_rm(all_of(list_env_predictors)) %>%
    recipes::update_role(-tidyselect::all_of(trait),-IDenv,-geno_ID, new_role = 'predictor') %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    #recipes::step_corr(recipes::all_predictors(),
    #                   skip = TRUE,
    #                   threshold = 0.95) %>%
    recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
  
  
  
  cat('Processing: recipe for the genomic-based kernel created!\n')
  
  ## KERNEL WITH INTERACTIONS BETWEEN SNPS AND ENVIRONMENTAL COVARIATES ##
  
  
  rec_GE <- recipes::recipe(~ . ,
                            data = training) %>%
    recipes::step_rm(-all_of(trait),-all_of(colnames(SNPs)),-all_of(list_env_predictors))  %>%
    recipes::step_rm('geno_ID') %>%
    recipes::step_interact(
      terms = ~ any_of(colnames(SNPs)):all_of(list_env_predictors),
      role = 'predictor'
    ) %>%
    recipes::step_rm(any_of(colnames(SNPs))) %>%
    recipes::step_rm(all_of(list_env_predictors)) %>%
    recipes::update_role(-all_of(trait), new_role = 'predictor') %>%
    recipes::update_role(all_of(trait), new_role = "outcome") %>%
    recipes::step_nzv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
  
  
  
  cat(
    paste(
      'Processing: recipe for the kernel with interactions between',
      'genomic and environmental predictors created!\n'
    )
  )
  
  split_processed <- structure(
    list(
      'training' = training,
      'test' = test,
      'rec_G' = rec_G,
      'rec_E' = rec_E,
      'rec_GE' = rec_GE
    ),
    class = 'stacking_reg_2'
  )
  
  
  
  return(split_processed)
  
  
}











#' @rdname stacking_reg_2
#' @aliases new_stacking_reg_2
#' @export
stacking_reg_2 <- function(split,
                           trait,
                           geno_data,
                           env_predictors,
                           info_environments,
                           geno_information,
                           use_selected_markers,
                           SNPs,
                           list_env_predictors,
                           include_env_predictors,
                           lat_lon_included,
                           year_included,
                           ...) {
  validate_stacking_reg_2(
    new_stacking_reg_2(
      split = split,
      trait = trait,
      geno_data = geno_data,
      env_predictors = env_predictors,
      info_environments = info_environments,
      geno_information = geno_information,
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


#' @rdname stacking_reg_2
#' @aliases new_stacking_reg_2
#' @export
validate_stacking_reg_2 <- function(x, ...) {
  trait <-
    as.character(x[['rec_G']]$term_info[which(x[['rec_G']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'stacking_reg_2')
  
  checkmate::assert_names(names(x),
                          must.include = c('training', 'test', 'rec_G', 'rec_E'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  
  return(x)
}
