#' Get train/test splits of the phenotypic MET dataset based on cv0.
#'
#' @description Get train/test splits of the phenotypic MET dataset based on a
#' number of random k-folds partitions determined by the user, according to the
#' type CV0. Creation of the train/test splits based on phenotypic data, so that
#' all the phenotypes from the same environment/year/site appear in the same 
#' fold, according to the type of the CV0 scheme.
#'
#' @param pheno_data \code{data.frame} Dataset containing phenotypic outcome
#'   data, as well as the predictor variables
#'
#' @param cv0_type \code{character} either 'leave-one-environment-out', 
#' 'leave-one-site-out', 'leave-one-year-out' or 'forward-prediction'
#'
#' @return a \code{list} which contains the train'test splits of the CV scheme.
#'   Each element of the list corresponds to a list with two elements:
#' \itemize{
#'   \item \code{data.frame} Dataset with all observations for the training set
#'   \item \code{data.frame} Dataset with all observations for the test set
#' }

predict_cv0 <-
  function(pheno_data,
           type = cv0_type) {
    
    if (cv0_type=='leave-one-year-out'){
      
      # Create data frame with unique names of year in the dataset
      
      unique_years <-
        as.character(unique(pheno_data[, 'year']))
      
      partition_data <- function(data, year) {
        
        training_data=data[data$year!=year,]
        test_data=data[data$year==year,]
        return(list(training_data, test_data))
        
      }
      
      train_test_splits <- map(
        unique_years,
        .f = function (x)
          partition_data(year=x, data = pheno_data)
      )
      
      return(train_test_splits)
    }

    if (cv0_type=='leave-one-environment-out'){
      
      # Create data frame with unique names of environments (YearxLoc) in the 
      # dataset
      
      unique_environments <-
        as.character(unique(pheno_data[, 'IDenv']))
      
      partition_data <- function(data, IDenv) {
        
        training_data=data[data$IDenv!=IDenv,]
        test_data=data[data$IDenv==IDenv,]
        return(list(training_data, test_data))
        
      }
      
      train_test_splits <- map(
        unique_environments,
        .f = function (x)
          partition_data(year=x, data = pheno_data)
      )
    
    return(train_test_splits)
    }
    
    if (cv0_type=='leave-one-site-out'){
      
      # Create data frame with unique names of location in the dataset
      
      unique_sites <-
        as.character(unique(pheno_data[, 'location']))
      
      partition_data <- function(data, location) {
        
        training_data=data[data$location!=location,]
        test_data=data[data$location==location,]
        return(list(training_data, test_data))
        
      }
      
      train_test_splits <- map(
        unique_environments,
        .f = function (x)
          partition_data(year=x, data = pheno_data)
      )
      
      return(train_test_splits)
    }
    
    
    if (cv0_type=='forward-prediction'){
      
      # Create data frame with unique names of year in the dataset
      
      pheno_data$year = as.numeric(as.character(pheno_data$year))
      unique_years <-
        as.numeric(unique(pheno_data[, 'year']))
      
      partition_data <- function(data, year) {
        
        training_data=data[data$year<year,]
        training_data$year = as.factor(training_data$year)
        test_data=data[data$year==year,]
        test_data$year = as.factor(test_data$year)
        return(list(training_data, test_data))
        
      }
      
      train_test_splits <- map(
        unique_years,
        .f = function (x)
          partition_data(year=x, data = pheno_data)
      )
      
      return(train_test_splits)
    
    
  }
}