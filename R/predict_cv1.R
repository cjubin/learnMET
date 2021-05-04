#' Get train/test splits of the phenotypic MET dataset based on CV1.
#' 
#' @description Get train/test splits of the phenotypic MET dataset based on a 
#' number of random k-folds partitions determined by the user, according to the
#' type CV1. Creation of the list of train/test splits based on phenotypic data,
#' so that all the phenotypes from the same line appear in same fold (prediction
#' of new lines never observed in any environment). 
#'
#' @param pheno_data \code{data.frame} Dataset containing phenotypic outcome
#'   data, as well as the predictor variables
#'
#' @param nb_folds \code{numeric} Number of folds in the CV process
#'
#' @param reps \code{numeric} Number of repeats of the k-folds CV
#'
#' @return a \code{list} which contains nb_folds x reps elements.
#'   Each element of the list corresponds to a list with two elements:
#' \itemize{
#'   \item \code{data.frame} Dataset with all observations for the training set
#'   \item \code{data.frame} Dataset with all observations for the test set
#' }


predict_cv1 <-
  function(pheno_data,
           nb_folds,
           reps,
           seed) {
    
    # Create data frame with unique names of lines
    
    unique_lines <-
      as.data.frame(unique(pheno_data[, 'geno_ID']))
    
    # Randomly assign lines to folds: k-fold cross-validation randomly splits
    # the lines into k folds of roughly equal size.
    # A resample of the analysis data consisted of K-1 of the folds while the
    # assessment set contains the final fold.
    
    set.seed(seed)
    
    lines_folds <-
      vfold_cv(data = unique_lines,
               v = nb_folds,
               repeats = reps)
    
    
    # Create the train/test splits of the phenotypic data so that all the
    # phenotypes from the same line appeared in same fold, according to the
    # resampling of lines previously done.
    
    partition_data <- function(splits, pheno) {
      
      training_lines <- analysis(splits)[, 1]
      test_lines <- assessment(splits)[, 1]
      
      training_data <- pheno[pheno$geno_ID %in% training_lines,]
      test_data <- pheno[pheno$geno_ID %in% test_lines,]
      return(list(training_data, test_data))
      
    }
    
    # Apply the function over the complete resampling object lines_folds (rset)
    
    
    train_test_splits <- map(
      lines_folds$splits,
      .f = function (x)
        partition_data(x, pheno = pheno_data)
    )
    
    return(train_test_splits)
    
  }
