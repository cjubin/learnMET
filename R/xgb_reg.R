#' @description
#' blabla
#' @title 
#' blabla
#' @name xgb_reg
#' @export
new_xgb_reg <- function(split,
                        trait,
                        geno_data,
                        env_predictors,
                        info_environments,
                        unique_EC_by_geno,
                        geno_information,
                        use_selected_markers,
                        SNPs,
                        list_env_predictors,
                        include_env_predictors,
                        lat_lon_included,
                        year_included) {
  if (class(split) != 'split') {
    stop('Class of x should be "split".')
  }
  
  if (class(split[['training']][, trait]) %in% c('integer')) {
    split[['training']][, trait] <-
      as.numeric(split[['training']][, trait])
    split[['test']][, trait] <- as.numeric(split[['test']][, trait])
    
  }
  
  
  ## GENOTYPIC DATA ##
  
  # Use of genotypic data: specified via parameter geno_information #
  
  if (geno_information == 'PCs') {
    # Processing of PCs: apply transformations calculated on the training set
    # on test set --> dimensionality reduction method
    
    cat('Processing: PCA transformation on the Training Set\n')
    pca_geno = apply_pca(split = split, geno = geno_data)
    
    # Merge in same data pheno and geno data for each train & test split
    
    training <-
      merge(split[[1]], pca_geno, by = 'geno_ID', all.x = T)
    
    test <-
      merge(split[[2]], pca_geno, by = 'geno_ID', all.x = T)
    
    cat('Processing: PCA transformation done\n')
    
  }
  
  
  # Add SNP covariates if they should be used
  
  if (use_selected_markers & geno_information != 'SNPs') {
    training <- merge(training, SNPs, by = 'geno_ID', all.x = T)
    test <- merge(test, SNPs, by = 'geno_ID', all.x = T)
    
  }
  
  
  ## ENVIRONMENTAL DATA ##
  # Add the environmental data
  
  # Two different cases:
  # Case 1: each environmental covariate is unique for a complete environment
  # (e.g. ECs are not computed individually for each variety within an
  # environment).
  
  if (include_env_predictors &
      !is.null(list_env_predictors) & !unique_EC_by_geno) {
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
  
  # Case 2: each environmental covariate is computed specifically for an
  # environment AND for a genotype (e.g. ECs are computed individually
  # for each variety within an environment).
  
  if (include_env_predictors &
      !is.null(list_env_predictors) & unique_EC_by_geno) {
    training <-
      merge(training,
            env_predictors[, c('IDenv', list_env_predictors)],
            by = c('IDenv', 'geno_ID'),
            all.x = T)
    test <-
      merge(test,
            env_predictors[, c('IDenv', list_env_predictors)],
            by = c('IDenv', 'geno_ID'),
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
      step_nzv(all_predictors(),-starts_with('PC')) %>%
      step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
    
    
    
    
    
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
  ), class = 'xgb_reg')
  
  
  
  
  return(split_processed)
  
  
}











#' @rdname xgb_reg
#' @aliases new_xgb_reg
#' @export
xgb_reg <- function(split,
                    trait,
                    geno_data,
                    env_predictors,
                    info_environments,
                    unique_EC_by_geno,
                    geno_information,
                    use_selected_markers,
                    SNPs,
                    list_env_predictors,
                    include_env_predictors,
                    lat_lon_included,
                    year_included) {
  validate_xgb_reg(
    new_xgb_reg(
      split = split,
      trait = trait,
      geno_data = geno_data,
      env_predictors = env_predictors,
      info_environments = info_environments,
      unique_EC_by_geno = unique_EC_by_geno,
      geno_information = geno_information,
      use_selected_markers = use_selected_markers,
      SNPs = SNPs,
      list_env_predictors = list_env_predictors,
      include_env_predictors = include_env_predictors,
      lat_lon_included = lat_lon_included,
      year_included = year_included
    )
  )
}


#' @rdname xgb_reg
#' @aliases new_xgb_reg
#' @export

validate_xgb_reg <- function(x) {
  trait <-
    as.character(x[['rec']]$term_info[which(x[['rec']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'xgb_reg')
  
  checkmate::assert_names(names(x), must.include = c('training', 'test', 'rec'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  checkmate::assert_class(x[['test']][, trait], 'numeric')
  
  
  return(x)
}
