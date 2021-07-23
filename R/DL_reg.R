#' Builds a recipe to process a split object (containing training 
#' and test sets) according to the configuration set by the user and assign it 
#' to a deep learning regression model for subsequent model fitting using a S3 
#' method dispatch.
#' 
#' @description
#' The function processes genomic information according to the option set by the
#' user. Predictor variables are standardized based on the training set, and
#' rhe use of recipes ensure that the same transformations are used on the test 
#' set.
#' based on the sets of predictors
#' geno info
#' @name new_DL_reg
#' deep learning reg
#' @name DL_reg
#' @export
new_DL_reg <- function(split,
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
  
  if (geno_information == 'PCs_G') {
    
    cat('Processing: PCs of the genomic relationship matrix\n')
    pcs_g_geno = apply_pcs_G(split = split, geno_data = geno_data,num_pcs=num_pcs)
    training = pcs_g_geno[[1]]
    test = pcs_g_geno[[2]]
    cat('Processing: PCs of the genomic relationship matrix done! \n')
    
  }
  
  if (geno_information == 'PCs') {
    
    cat('Processing: PCA transformation on the scaled marker dataset\n')
    pca_geno = apply_pca(split = split, geno_data = geno_data,num_pcs=num_pcs)
    training = pca_geno[[1]]
    test = pca_geno[[2]]
    cat('Processing: PCA transformation done\n')
    
  }
  
  if (geno_information == 'SNPs') {
    
    geno_data$geno_ID = row.names(geno_data)
    
    geno_training = geno_data[geno_data$geno_ID%in%unique(split[[1]][,'geno_ID']),]
    geno_training = unique(geno_training)
    geno_test =  geno_data[geno_data$geno_ID%in%unique(split[[2]][,'geno_ID']),]
    geno_test = unique(geno_test)
    
    rec_snps <- recipe(~ . ,
                       data = geno_training) %>%
      update_role(geno_ID, new_role = 'outcome') %>%
      step_normalize(all_numeric_predictors()) 
    
    rec_snps <- prep(rec_snps,training = geno_training,strings_as_factors = FALSE)
    
    snps_data_tr <- bake(rec_snps,new_data = geno_training)
    snps_data_te <- bake(rec_snps, new_data = geno_test)
    
    training <-
      merge(split[[1]], snps_data, by = 'geno_ID', all.x = T)
    
    test <-
      merge(split[[2]], snps_data, by = 'geno_ID', all.x = T)
    
    
  }
  
  # Add SNP covariates if they should be used
  
  if (use_selected_markers & geno_information != 'SNPs') {
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
  ), class = 'DL_reg')
  
  
  
  
  return(split_processed)
  
  
}











#' @rdname DL_reg
#' @aliases new_DL_reg
#' @export
DL_reg <- function(split,
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
  validate_DL_reg(
    new_DL_reg(
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


#' @rdname DL_reg
#' @aliases new_DL_reg
#' @export

validate_DL_reg <- function(x,...) {
  trait <-
    as.character(x[['rec']]$term_info[which(x[['rec']]$term_info[, 3] == 'outcome'), 'variable'])
  
  checkmate::assert_class(x, 'DL_reg')
  
  checkmate::assert_names(names(x), must.include = c('training', 'test', 'rec'))
  
  checkmate::assert_class(x[['training']][, trait], 'numeric')
  
 
  
  return(x)
}
