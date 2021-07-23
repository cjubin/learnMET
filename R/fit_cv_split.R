#' S3 method used to fit an object of class `xgb_reg`, `xgb_ordinal`,
#' `svm_stacking_reg` or `DL_reg`.
#' 
#' @description
#' S3 dispatching method for objects of class `xgb_reg`, `xgb_ordinal`,
#' `svm_stacking_reg` or `DL_reg`.
#' 
#' @name fit_cv_split
#'
#' @param object an object of class `xgb_reg`, `xgb_ordinal`,
#' `svm_stacking_reg` or `DL_reg`.
#' 
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
fit_cv_split <- function(object, ...) {
  UseMethod("fit_cv_split")
}


#' @rdname fit_cv_split
#' @export
fit_cv_split.default <- function(x, ...) {
  stop('not implemented')
  
}
