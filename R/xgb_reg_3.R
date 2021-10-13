#' Processing of a split object to get data ready to be used and fitted with
#' a `xgb_reg_3` (gradient boosted tree) regression model.
#'
#' @description
#' The function processes a split object (training + test sets), according to
#' the configuration set by the user. For instance, genomic information is 
#' incorporated according to the option set by the user.  A list of specific
#' environmental covariables to use can be provided.\cr
#' 
#' A recipe is created using the package `recipes`, to specify additional 
#' preprocessing steps, such as standardization based on the training set, with
#' same transformations used on the test set. Variables with null variance are
#' removed. If year effect is included, it is converted to dummy variables.\cr
#' Further fitting on the training set with a gradient boosting model (see 
#' function [fit_cv_split.xgb_reg_3()])).
#' 
#' \strong{This prediction method can be very slow according to the number of SNPs variables used!}
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
#' @return A `list` object of class `xgb_reg_3` with the following items:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set after partial processing}
#'     \item{test}{\code{data.frame} Test set after partial processing}
#'     \item{rec}{A \code{recipe} object, specifying the remaining processing
#'     steps which are implemented when a model is fitted on the training set
#'     with a recipe.}
#'   }
#'     
#' @references
#' \insertRef{wickham2019welcome}{learnMET}
#' \insertRef{tidymodels}{learnMET}
#'
#' @name xgb_reg_3
#' @export
new_xgb_reg_3 <- function(split = NULL,
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
  
  
  ## GENOTYPIC DATA ##
  
  # Use of genotypic data: use of all SNPs as predictor variables ! CAN BE VERY SLOW ! #
  
    geno$geno_ID = row.names(geno)
    
    geno_training = geno[geno$geno_ID%in%unique(split[[1]][,'geno_ID']),]
    geno_training = unique(geno_training)
    geno_test =  geno[geno$geno_ID%in%unique(split[[2]][,'geno_ID']),]
    geno_test = unique(geno_test)
     
    rec_snps <- recipe(~ . ,
                  data = geno_training) %>%
      update_role(geno_ID, new_role = 'outcome') %>%
      step_normalize(all_numeric_predictors()) 
    
    rec_snps <- recipes::prep(rec_snps,training = geno_training,strings_as_factors = FALSE)
    
    snps_data_tr <- recipes::bake(rec_snps,new_data = geno_training)
    snps_data_te <- recipes::bake(rec_snps, new_data = geno_test)
    
    training <-
      merge(split[[1]], snps_data_tr, by = 'geno_ID', all.x = T)
    
    test <-
      merge(split[[2]], snps_data_te, by = 'geno_ID', all.x = T)
   
    
  
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
    
    rec <- recipe(~ . ,
                  data = training) %>%
      update_role(trait, new_role = 'outcome') %>%
      update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      step_rm(location) %>%
      step_rm(geno_ID) %>%
      update_role(-trait,-IDenv, new_role = 'predictor') %>%
      step_dummy(year, preserve = F, one_hot = TRUE) %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric(),-all_outcomes())
    
    
    
    
    
  }
  else if (!lat_lon_included &
           year_included &
           length(unique(as.character(training$year))) > 1) {
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipe(~ . ,
                  data = training) %>%
      update_role(trait, new_role = 'outcome') %>%
      update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      step_rm(location) %>%
      step_rm(geno_ID) %>%
      update_role(-trait,-IDenv, new_role = 'predictor') %>%
      step_dummy(year, preserve = F, one_hot = TRUE) %>%
      step_nzv(all_predictors(),-starts_with('PC')) %>%
      step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
    
    
    
    
  }
  
  else if ((lat_lon_included &
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
    
    rec <- recipe(~ . ,
                  data = training) %>%
      update_role(trait, new_role = 'outcome') %>%
      update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      step_rm(location) %>%
      step_rm(geno_ID) %>%
      step_rm(year) %>%
      update_role(-trait,-IDenv, new_role = 'predictor') %>%
      step_nzv(all_predictors(),-starts_with('PC')) %>%
      step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
    
    
    
    
  }
  else{
    # Create recipe to define the processing of the training & test set.
    
    rec <- recipe(~ . ,
                  data = training) %>%
      update_role(trait, new_role = 'outcome') %>%
      update_role(IDenv, location, geno_ID, new_role = "id variable") %>%
      step_rm(location) %>%
      step_rm(geno_ID) %>%
      step_rm(year) %>%
      update_role(-trait,-IDenv, new_role = 'predictor') %>%
      step_nzv(all_predictors(),-starts_with('PC')) %>%
      step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
    
    
    
    
  }
  cat(
    'Incorporating selected predictors & Data processing for one train/test split of the CV scheme: Done!\n'
  )
  
  split_processed <- structure(list(
    "training" = training,
    "test" = test,
    "rec" = rec
  ), class = 'xgb_reg_3')
  
  
  
  
  return(split_processed)
  
  
}






#' @rdname xgb_reg
#' @aliases new_xgb_reg
#' @export
xgb_reg <- function(split,
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
  
  validate_xgb_reg(
    new_xgb_reg(
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


#' @rdname xgb_reg
#' @aliases new_xgb_reg
#' @export

validate_xgb_reg <- function(x,...) {
  trait <-
    as.character(x[['rec']]$term_info[which(x[['rec']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'xgb_reg_3')
  
  checkmate::assert_names(names(x), must.include = c('training', 'test', 'rec'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  
  
  return(x)
}
