


predict_cv1 <-
  function(pheno_all_data,
           nb_folds = nb_folds_cv1,
           reps = repeats_cv1) {
    # Create data frame with unique names of lines
    
    unique_lines <-
      as.data.frame(unique(pheno_all_data[, 'geno_ID']))
    
    # Randomly assign lines to folds: k-fold cross-validation randomly splits
    # the lines into k folds of roughly equal size.
    # A resample of the analysis data consisted of K-1 of the folds while the
    # assessment set contains the final fold.
    
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
      
      training_data <- pheno[pheno$geno_ID %in% training_lines,-c('geno_ID')]
      test_data <- pheno[pheno$geno_ID %in% test_lines,,-c('geno_ID')]
      return(list(training_data, test_data))
      
    }
    
    # Apply the function over the complete resampling object lines_folds (rset)
    
    
    train_test_splits <- map(
      lines_folds$splits,
      .f = function (x)
        partition_data(x, pheno = pheno_all_data)
    )
    
    
  }
