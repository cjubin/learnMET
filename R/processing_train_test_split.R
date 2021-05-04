#' Processing and selecting predictors to fit the model
#'
#' @param split
#' @param list_env_predictors
#' @param env_predictors
#' @param lat_lon_included
#' @param year_included
#' @param use_selected_markers
#'
#' @return 
#'
#'
#'


processing_train_test_split <-
  function(split,
           geno_data,
           geno_information,
           use_selected_markers,
           SNPs,
           list_env_predictors,
           include_env_predictors,
           lat_lon_included,
           year_included) {
    
    ## GENOTYPIC DATA ## 
    
    # Use of genotypic data: specified via parameter geno_information #
    
    if (geno_information== 'PCs') {
      # Processing of PCs: apply transformations calculated on the training set
      # on test set --> dimensionality reduction method
      
      cat('Processing: PCA transformation on the Training Set\n')
      pca_geno = apply_pca(split = split,geno=geno_data)
      
      # Merge in same data pheno and geno data for each train & test split
      
      training <-
        merge(split[[1]], pca_geno, by = 'geno_ID', all.x = T)
      
      test <-
        merge(split[[2]], pca_geno, by = 'geno_ID', all.x = T)
      
      cat('Processing: PCA transformation done\n')
      
    }
    
   
    # Add SNP covariates if they should be used
    
    if (use_selected_markers) {
      
      training <- merge(training, SNPs, by = 'geno_ID', all.x = T)
      
      test <- merge(test, SNPs, by = 'geno_ID', all.x = T)
      
    }
    
    
    ## ENVIRONMENTAL DATA ## 
    
    if (include_env_predictors & !is.null(list_env_predictors)) {
      
      # Add the environmental data
      
      training <-
        merge(training,
              METData$env_data[, c('IDenv', list_env_predictors)],
              by = 'IDenv',
              all.x = T)
      test <-
        merge(test,
              METData$env_data[, c('IDenv', list_env_predictors)],
              by = 'IDenv',
              all.x = T)
      
    }
    
    if (lat_lon_included & year_included & length(unique(as.character(training$year)))>1) {
      # Add longitude/latitude data for each train & test split
      
      training <-
        merge(training,
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
              by = 'IDenv',
              all.x = T)
      test <-
        merge(test,
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
              by = 'IDenv',
              all.x = T)

      
      
      # Create recipe to define the processing of the training & test set.
      
      rec <- recipe(~ . ,
                    data = training) %>%
        update_role(trait, new_role = 'outcome') %>%
        update_role(IDenv, new_role = "id variable") %>%
        step_rm(location) %>%
        step_rm(geno_ID) %>%
        update_role(-trait,-IDenv, new_role = 'predictor') %>%
        step_dummy(year, preserve = F, one_hot = TRUE) %>%
        step_nzv(all_predictors(),-starts_with('PC'))%>%
        step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) 
        
      
     
      
      
    }
    else if (!lat_lon_included &
             year_included & length(unique(as.character(training$year)))>1) {
      
      
      # Create recipe to define the processing of the training & test set.
      
      rec <- recipe(~ . ,
                    data = training) %>%
        update_role(trait, new_role = 'outcome') %>%
        update_role(IDenv, new_role = "id variable") %>%
        step_rm(location) %>%
        step_rm(geno_ID) %>%
        update_role(-trait,-IDenv, new_role = 'predictor') %>%
        step_dummy(year, preserve = F, one_hot = TRUE) %>%
        step_nzv(all_predictors(),-starts_with('PC')) %>%
        step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC')) 
      
     
      
      
    }
    
    else if ((lat_lon_included &
             !year_included) | (lat_lon_included &
                                length(unique(as.character(training$year)))<2) ) {
      # Add longitude/latitude data for each train & test split
      
      training <-
        merge(training,
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
              by = 'IDenv',
              all.x = T)
      test <-
        merge(test,
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
              by = 'IDenv',
              all.x = T)
      
      
      # Create recipe to define the processing of the training & test set.
      
      rec <- recipe(~ . ,
                    data = training) %>%
        update_role(trait, new_role = 'outcome') %>%
        update_role(IDenv, new_role = "id variable") %>%
        step_rm(location) %>%
        step_rm(geno_ID) %>%
        step_rm('year') %>%
        update_role(-trait,-IDenv, new_role = 'predictor') %>%
        step_nzv(all_predictors(),-starts_with('PC')) %>%
        step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
      
     
      
      
    }
    else{
      
      # Create recipe to define the processing of the training & test set.
      
      rec <- recipe(~ . ,
                    data = training) %>%
        update_role(trait, new_role = 'outcome') %>%
        update_role(IDenv, new_role = "id variable") %>%
        step_rm(location) %>%
        step_rm(geno_ID) %>%
        step_rm('year') %>%
        update_role(-trait,-IDenv, new_role = 'predictor') %>%
        step_nzv(all_predictors(),-starts_with('PC')) %>%
        step_normalize(all_numeric(),-all_outcomes(), -starts_with('PC'))
      
     
      
      
    }
    cat('Incorporating selected predictors & Data processing for one train/test split of the CV scheme: Done!\n')
    
    return(list(training, test, rec))
    
  }
