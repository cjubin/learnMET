#' Processing and selecting predictors to fit the model
#'
#' @param split
#' @param list_env_predictors
#' @param env_predictors
#' @param lat_lon_included
#' @param year_included
#' @param use_selected_markers
#' @param geno threshold 0.95 correlation
#'
#' @return
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'


processing_train_test_split_kernel <-
  function(split,
           geno_data,
           geno_information,
           use_selected_markers,
           SNPs,
           list_env_predictors,
           include_env_predictors,
           lat_lon_included = F,
           year_included = F,
           ...) {
    
    
    geno_data$geno_ID = row.names(geno_data)
    
    # Merge in same data.frame pheno and geno data for each train & test split
    
    training <-
      merge(split[[1]], geno_data, by = 'geno_ID', all.x = T)
    
    test <-
      merge(split[[2]], geno_data, by = 'geno_ID', all.x = T)
    
    # Add the environmental data
    
    if (include_env_predictors & !is.null(list_env_predictors)) {
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
    

    
    ## ENVIRONMENTAL-BASED KERNEL ##
    
    if (lat_lon_included &
        year_included &
        length(unique(as.character(training$year))) > 1) {
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
      
      rec_E <- recipes::recipe(~ . ,
                               data = training) %>%
        recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
        recipes::update_role(IDenv, new_role = "id variable") %>%
        recipes::step_rm(all_of(colnames(geno_data))) %>%
        recipes::step_rm(location) %>%
        recipes::step_rm(colnames(geno_data)) %>%
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
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
              by = 'IDenv',
              all.x = T)
      test <-
        merge(test,
              METData$info_environments[, c('IDenv', 'longitude', 'latitude')],
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
      recipes::update_role(tidyselect::all_of(trait), new_role = 'outcome') %>%
      recipes::update_role(IDenv, new_role = "id variable") %>%
      recipes::step_interact(
        terms = ~ all_of(colnames(SNPs)[colnames(SNPs) %notin% 'geno_ID']):all_of(list_env_predictors),
        role = 'predictor')  %>%
      recipes::step_rm(all_of(colnames(geno_data))) %>%
      recipes::step_rm(all_of(list_env_predictors)) %>%
      recipes::step_rm(location) %>%
      recipes::step_rm(year) %>%
      recipes::update_role(-all_of(trait),-IDenv, new_role = 'predictor') %>%
      recipes::step_nzv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric(),-recipes::all_outcomes())
    
    
    cat(paste('Processing: recipe for the kernel with interactions between',
              'genomic and environmental predictors created!\n'))
    
    return(
      list(
        'training' = training,
        'test' = test,
        'rec_G' = rec_G,
        'rec_E' = rec_E,
        'rec_GE' = rec_GE
      )
    )
    
    
    
  }
