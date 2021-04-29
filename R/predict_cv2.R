#' Get train/test splits of the phenotypic MET dataset based on CV2.
#' 
#' @description Get train/test splits of the phenotypic MET dataset based on a 
#' number of random k-folds partitions determined by the user, according to the
#' type CV2. Creation of the train/test splits based on phenotypic data, so that
#' all the YearxLocation phenotypic observations from the phenotypic MET dataset
#' are assigned randomly to k-fold partitions (prediction of incomplete field 
#' trials).
#'
#' @param pheno_data \code{data.frame} Dataset containing phenotypic outcome
#'   data, as well as the predictor variables
#'
#' @param nb_folds \code{numeric} Number of folds in the CV process
#'
#' @param reps \code{numeric} Number of repeats of the k-folds CV
#'
#' @return a \code{list} which contains nb_folds x reps elements
#'   Each element of the list corresponds to a list with two elements:
#' \itemize{
#'   \item \code{data.frame} Dataset with all observations for the training set
#'   \item \code{data.frame} Dataset with all observations for the test set
#' }

predict_cv2 <-
  function(pheno_data,
           nb_folds = nb_folds_cv2,
           reps = repeats_cv2) {
     
    # Randomly assign phenotypic observations to folds: k-fold cross-validation 
    # randomly splits the lines into k folds of roughly equal size.
    # A resample of the analysis data consisted of K-1 of the folds while the
    # assessment set contains the final fold.
    
    lines_folds <-
      vfold_cv(data = pheno_data,
               v = nb_folds,
               repeats = reps)
    
    
    partition_data <- function(splits, pheno) {
      
      training_data <- analysis(splits)
      test_data <- assessment(splits)
      
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
