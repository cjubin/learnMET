#' @description
#' blabla
#' @title ge
#' @name svm_stacking_reg
#' @export
new_svm_stacking_reg <- function(split,
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
                                 year_included,
                                 ...) {
  
  if (class(split) != 'split') {
    stop('Class of x should be "split".')
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
      merge(
        training,
        env_predictors[, c('IDenv', list_env_predictors)],
        by = c('IDenv', 'geno_ID'),
        all.x = T
      )
    test <-
      merge(test,
            env_predictors[, c('IDenv', list_env_predictors)],
            by = c('IDenv', 'geno_ID'),
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
  
  if (use_selected_markers) {
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
      class = 'svm_stacking_reg'
    )
    
  } else {
    split_processed <- structure(list(
      'training' = training,
      'test' = test,
      'rec_G' = rec_G,
      'rec_E' = rec_E
    ),
    class = 'svm_stacking_reg')
  }
  
  
  return(split_processed)
  
  
}











#' @rdname svm_stacking_reg
#' @aliases new_svm_stacking_reg
#' @export
svm_stacking_reg <- function(split,
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
                             year_included,
                             ...) {
  validate_svm_stacking_reg(
    new_svm_stacking_reg(
      split=split,
      trait=trait,
      geno_data=geno_data,
      env_predictors = env_predictors,
      info_environments = info_environments,
      unique_EC_by_geno = unique_EC_by_geno,
      geno_information=geno_information,
      use_selected_markers=use_selected_markers,
      SNPs=SNPs,
      list_env_predictors=list_env_predictors,
      include_env_predictors=include_env_predictors,
      lat_lon_included=lat_lon_included,
      year_included=year_included,
      ...
    )
  )
}


#' @rdname svm_stacking_reg
#' @aliases new_svm_stacking_reg
#' @export
validate_svm_stacking_reg <- function(x,...) {
  trait <-
    as.character(x[['rec_G']]$term_info[which(x[['rec_G']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'svm_stacking_reg')
  
  checkmate::assert_names(names(x),
                          must.include = c('training', 'test', 'rec_G', 'rec_E'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
  
  return(x)
}
