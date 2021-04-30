#' Processing and selecting predictors to fit the model
#'
#' @param
#'
#'
#'
#'


processing_train_test_split <-
  function(split,
           list_env_predictors,
           include_env_predictors,
           lat_lon_included,
           year_included,
           use_selected_markers) {
    # Merge in same data pheno and geno data for each train & test split
    
    training <-
      merge(split[[1]], split[[3]], by = 'geno_ID', all.x = T)
    
    test <-
      merge(split[[2]], split[[3]], by = 'geno_ID', all.x = T)
    
    # Add SNP covariates if they should be used
    
    if (use_selected_markers) {
      training <- merge(training, SNPs, by = 'geno_ID', all.x = T)
      
      test <- merge(test, SNPs, by = 'geno_ID', all.x = T)
      
    }
    
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
    cat('Incorporating selected predictors & Data processing for one train/test split of the CV scheme: Done!')
    
    return(list(training, test, rec))
    
  }
