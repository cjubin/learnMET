new_xgb_multiclass_factor <- function(processed_splits, trait) {
  
  
  if (class(processed_splits[[1]]) != 'split_processed') {
    stop('Class of x should be "split_processed".')
  }
  
  if (class(processed_splits[[1]][['training']][, trait]) %in% c('integer')) {
    to_factor <-
      function(y) {
        y[['training']][, trait] <- as.factor(y[['training']][, trait])
        y[['test']][, trait] <- as.factor(y[['test']][, trait])
        return(y)
      }
    processed_splits_transformed <-
      lapply(processed_splits, function(x) {
        to_factor(y = x)
      })
    
  }
  
  xgb_multiclass_factor <-
    structure(processed_splits_transformed, class = 'xgb_multiclass_factor')
  
  
  return(xgb_multiclass_factor)
  
  
}


xgb_multiclass_factor <- function(processed_splits, trait) {
  validate_xgb_multiclass_factor(new_xgb_multiclass_factor(processed_splits, trait))
}


validate_xgb_multiclass_factor <- function(x) {
  checkmate::assert_class(x, 'xgb_multiclass_factor')
  
  trait <- as.character(x[[1]][['rec']]$term_info[which(x[[1]][['rec']]$term_info[,3]=='outcome'),'variable'])
  
  check_subelements_names <-
    function(y) {
      checkmate::test_names(names(y),
                            must.include = c('training', 'test', 'rec'))
    }
  
  all_splits_check_names <- unlist(lapply(x, function(x) {
    check_subelements_names(y=x)
  }))
  
  if (all(all_splits_check_names)!=TRUE){stop('Wrong names for each subelement of the xgb_multiclass_factor class object.')}
  
  check_new_factor_class <-
    function(y) {
      checkmate::test_class(y[['training']][,trait],
                            'factor')
    }
  
  all_splits_check_names <- unlist(lapply(x, function(x) {
    check_new_factor_class(y=x)
  }))
  
  if (all(all_splits_check_names)!=TRUE){stop('Trait has not been correctly converted to factor in the xgb_multiclass_factor class object.')}
  
  
  return(x)
}

