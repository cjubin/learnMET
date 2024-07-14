#' Get train/test splits of the phenotypic MET dataset based on CV0.
#'
#' @description Get train/test splits of the phenotypic MET dataset based on a
#' number of random k-folds partitions determined by the user, according to the
#' type CV00. Creation of the list of train/test splits based on phenotypic
#' data, so that all the phenotypes from the same environment/year/site appear
#' in the same fold, according to the type of the CV00 scheme. In addition to
#' CV0 scheme, information on lines present in the test set evaluated in other
#' environments are removed from the training set --> prediction of new
#' genotypes in new environments.
#'
#' @param pheno_data \code{data.frame} Dataset containing phenotypic outcome
#'   data, as well as the predictor variables.
#'
#' @param cv0_type \code{character} either `leave-one-environment-out`,
#'   `leave-one-site-out`, `leave-one-year-out` or `forward-prediction`.
#'
#' @param nb_folds \code{integer}
#'
#' @return a \code{cv_object} object which contains the train/test splits of the
#'   CV scheme. Each element of the object corresponds to a `split` object with
#'   two elements:
#'   \describe{
#'     \item{training}{\code{data.frame} Dataset with all observations for the
#'      training set.}
#'     \item{test}{\code{data.frame} Dataset with all observations for the test
#'      set.}
#'   }
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @references
#' \insertRef{jarquin2017increasing}{learnMET}
#' \insertRef{jarquin2014reaction}{learnMET}
#' @export

predict_cv00_5folds <-
  function(pheno_data,
           cv0_type,
           nb_folds_cv00) {
    if (cv0_type == 'leave-one-year-out') {
      # Create data frame with unique names of year in the dataset
      
      unique_years <-
        as.character(unique(pheno_data[, 'year']))
      
      partition_data <- function(data, year) {
        test_data = data[data$year == year,]
        training_data = data[data$year != year,]
        lines_test_set = unique(test_data$geno_ID)
        training_data = training_data[training_data$geno_ID %notin% lines_test_set,]
        
        split <-
          list("training_data" = training_data, "test_data" = test_data)
        class(split) <- c('split')
        names(split) <- c('training', 'test')
        return(split)
        
      }
      
      train_test_splits <- purrr::map(
        unique_years,
        .f = function (x)
          partition_data(year = x, data = pheno_data)
      )
      class(train_test_splits) <- c('cv_object')
      return(train_test_splits)
    }
    
    if (cv0_type == 'leave-one-environment-out') {
      # Create data frame with unique names of environments (YearxLoc) in the
      # dataset
      unique_environments <-
        as.character(unique(pheno_data[, 'IDenv']))
      
      
      partition_data_per_environment <-
        function(data, IDenv, nb_folds) {
          table_genotypes <-
            as.data.frame(unique(data[which(data$IDenv %in% IDenv), ]$geno_ID))
          colnames(table_genotypes)[1] <- 'geno_ID'
          set.seed(100)
          folds <-
            rsample::vfold_cv(table_genotypes, v = nb_folds, repeats = 1)
          
          
          within_environment_fold_split_cv00  <-
            function(nb_fold,
                     table_genotypes,
                     data,
                     IDenv) {
              training_data = data[data$IDenv != IDenv,]
              test_data = data[data$IDenv == IDenv,]
              
              sample_genotypes_training_set <-
                table_genotypes$geno_ID[folds$splits[[nb_fold]][[2]]]
              sample_genotypes_test_set <-
                table_genotypes$geno_ID[table_genotypes$geno_ID %notin% sample_genotypes_training_set]
              
              training_data <-
                training_data[which(training_data$geno_ID %in% sample_genotypes_training_set), ]
              test_data <-
                test_data[which(test_data$geno_ID %in% sample_genotypes_test_set), ]
              
              
              split <-
                list("training_data" = training_data, "test_data" = test_data)
              class(split) <- c('split')
              names(split) <- c('training', 'test')
              return(split)
            }
          
          train_test_splits_within_env <- purrr::map(
            .x = c(1:5),
            .f = function (x)
              within_environment_fold_split_cv00(
                nb_fold = x,
                data = pheno_data,
                table_genotypes = table_genotypes,
                IDenv = IDenv
              )
          )
          
          
          
        }
      
      train_test_splits <- purrr::map(
        .x = unique_environments,
        .f = function (x)
          partition_data_per_environment(IDenv = x, data = pheno_data, nb_folds = nb_folds_cv00))
      flatten_list <- list()
      n_tr_te_splits <- nb_folds_cv00 * length(unique_environments)
      n = 1
      for (env in 1:length(unique_environments)){
        for (f in 1:nb_folds_cv00){
          flatten_list[[n]] <- train_test_splits[[env]][[f]]
          n <- n+1
        }
      }
      class(flatten_list) <- c('cv_object')
      return(flatten_list)
    }
    
    if (cv0_type == 'leave-one-site-out') {
      # Create data frame with unique names of location in the dataset
      
      unique_sites <-
        as.character(unique(pheno_data[, 'location']))
      
      partition_data <- function(data, location) {
        training_data = data[data$location != location,]
        test_data = data[data$location == location,]
        lines_test_set = unique(test_data$geno_ID)
        training_data = training_data[training_data$geno_ID %notin% lines_test_set,]
        
        split <-
          list("training_data" = training_data, "test_data" = test_data)
        class(split) <- c('split')
        names(split) <- c('training', 'test')
        return(split)
        
      }
      
      train_test_splits <- map(
        unique_environments,
        .f = function (x)
          partition_data(location = x, data = pheno_data)
      )
      class(train_test_splits) <- c('cv_object')
      return(train_test_splits)
    }
    
    
    if (cv0_type == 'forward-prediction') {
      # Create data frame with unique names of year in the dataset
      
      unique_years = unique(as.numeric(as.character(pheno_data$year)))
      
      unique_years <- unique_years[-which.min(unique_years)]
      
      pheno_data$year <- as.numeric(as.character(pheno_data$year))
      
      partition_data <- function(data, year) {
        training_data = data[data$year < year,]
        training_data$year = as.factor(training_data$year)
        test_data = data[data$year == year,]
        test_data$year = as.factor(test_data$year)
        lines_test_set = unique(test_data$geno_ID)
        training_data = training_data[training_data$geno_ID %notin% lines_test_set,]
        
        split <-
          list("training_data" = training_data, "test_data" = test_data)
        class(split) <- c('split')
        names(split) <- c('training', 'test')
        return(split)
        
      }
      
      train_test_splits <- map(
        unique_years,
        .f = function (x)
          partition_data(year = x, data = pheno_data)
      )
      
      train_test_splits <-
        train_test_splits[lengths(train_test_splits) != 0]
      class(train_test_splits) <- c('cv_object')
      return(train_test_splits)
      
      
    }
  }