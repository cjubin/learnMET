

predict_cv2 <-
  function(pheno_all_data,
           nb_folds = nb_folds_cv2,
           reps = repeats_cv2) {
     
    # Randomly assign phenotypic observations to folds: k-fold cross-validation randomly splits
    # the lines into k folds of roughly equal size.
    # A resample of the analysis data consisted of K-1 of the folds while the
    # assessment set contains the final fold.
    
    lines_folds <-
      vfold_cv(data = pheno_all_data,
               v = nb_folds,
               repeats = reps)
    
    # Create the train/test splits of the phenotypic data so that all the
    # phenotypes from the same line appeared in same fold, according to the
    # resampling of lines previously done.
    
    partition_data <- function(splits, pheno) {
      
      training_data <- analysis(splits)
      test_data <- assessment(splits)
      
      return(list(training_data, test_data))
      
    }
    
    # Apply the function over the complete resampling object lines_folds (rset)
    
    
    train_test_splits <- map(
      lines_folds$splits,
      .f = function (x)
        partition_data(x, pheno = pheno_all_data)
    )
    
    
  }
