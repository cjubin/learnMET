#' S3 method used to fit an object of class  `rf_reg_1`, `rf_reg_2`, 
#' `rf_reg_3`, `xgb_reg_1`, `xgb_reg_2`, `xgb_reg_3`,`DL_reg`,`DL_reg_1`,
#' `DL_reg_2`,`DL_reg_3`,`stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' 
#' @description
#' S3 dispatching method for objects of class `rf_reg_1`, `rf_reg_2`, 
#' `rf_reg_3`, `xgb_reg_1`, `xgb_reg_2`, `xgb_reg_3`,`DL_reg`,`DL_reg_1`,
#' `DL_reg_2`,`DL_reg_3`,`stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' 
#' @name fit_split
#'
#' @param object an object of class `rf_reg_1`, `rf_reg_2`, `rf_reg_3`, 
#' `xgb_reg_1`, `xgb_reg_2`, `xgb_reg_3`,`DL_reg`,`DL_reg_1`,`DL_reg_2`,
#' `DL_reg_3`,`stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @export
fit_split <- function(object, ...) {
  UseMethod("fit_split")
}


#' @rdname fit_split
#' @export
fit_split.default <- function(object, ...) {
  stop('not implemented')
  
}
