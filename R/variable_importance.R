#' S3 method used to fit an object of class `xgb_reg`, `DL_reg`,
#' `stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' @description
#' S3 dispatching method for objects of class `xgb_reg`, `DL_reg`,
#' `stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' @name variable_importance_split
#'
#' @param object an object of class `xgb_reg`, `DL_reg`,
#' `stacking_reg_1`, `stacking_reg_2` or `stacking_reg_3`.
#' 
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
variable_importance_split <- function(object, ...) {
  UseMethod("variable_importance_split")
}


#' @rdname variable_importance_split
#' @export
variable_importance_split.default <- function(x, ...) {
  stop('not implemented')
  
}

#' @rdname variable_importance_split
#' @export
variable_importance_split.DL_reg <- function(fitted_obj_for_vip) {
   
  model <- fitted_obj_for_vip$model
  trait <- fitted_obj_for_vip$trait
  y_train <- fitted_obj_for_vip$y_train
  x_train <- fitted_obj_for_vip$x_train
  
  # create custom predict function
  pred_wrapper <- function(model, newdata)  {
    results <- predict(model, x = newdata) %>%
       as.vector()
    return(results)
    
  }
  print(pred_wrapper(model=model,newdata=x_train))
  
  explainer <- DALEX::explain(
    model = model,
    data = x_train,
    y = as.numeric(y_train),
    predict_function = pred_wrapper,
    type = 'regression',
    label = "DL_model_vip"
  ) 
  ranking_vip <- model_parts(explainer, N=NULL, loss_function = loss_root_mean_square,B=50)
  
  ranking_vip <- ranking_vip[,c('variable','mean_dropout_loss')]
  colnames(ranking_vip) <- c('Variable','Importance')
  
  return(ranking_vip)
}

#' @rdname variable_importance_split
#' @export
variable_importance_split.xgb_reg <- function(fitted_obj_for_vip) {
  
  
  # Obtain the variable importance with the gain metric
  
  model <- fitted_obj_for_vip$model
  trait <- fitted_obj_for_vip$trait
  y_train <- fitted_obj_for_vip$y_train
  x_train <- fitted_obj_for_vip$x_train
  
  predictors <- model %>%
    fit(data = x_train) %>%
    pull_workflow_fit()
  
  predictors <- predictors$fit$feature_names
  
  
  variable_importance_vip <- model %>%
    fit(data = x_train) %>%
    pull_workflow_fit() %>%
    vip::vip(num_features = length(predictors))
  
  
  ranking_vip <- as.data.frame(variable_importance_vip$data)
  
  if (length(predictors[which(predictors %notin% ranking_vip$Variable)])
      >
      0) {
    remaining <-
      cbind(as.vector(predictors[which(predictors %notin% ranking_vip$Variable)]), as.numeric(0))
    colnames(remaining) <- colnames(ranking_vip)
    ranking_vip <- rbind(ranking_vip, remaining)
  }
  ranking_vip$Importance <- as.numeric(ranking_vip$Importance)
  
  
  
  return(ranking_vip)
}

#' @rdname variable_importance_split
#' @export
variable_importance_split.stacking_reg_1 <- function(fitted_obj_for_vip) {
  
  model <- fitted_obj_for_vip$model
  trait <- fitted_obj_for_vip$trait
  y_train <- as.numeric(fitted_obj_for_vip$y_train)
  x_train <- fitted_obj_for_vip$x_train

  
  # create custom predict function
  pred_wrapper <- function(model, newdata)  {
    results <- as.data.frame(model %>% predict(new_data = newdata))
    
    return(results)
  }

  
  explainer <- DALEXtra::explain_tidymodels(
    model = model,
    data = x_valid,
    y = y_valid,
    predict_function = pred_wrapper,
    label = "stacked_model_vip"
  ) 
  ranking_vip <- model_parts(explainer, N=NULL, loss_function = loss_root_mean_square,B=50)
  
  ranking_vip <- ranking_vip[,c('variable','mean_dropout_loss')]
  colnames(ranking_vip) <- c('Variable','Importance')
  
  
  return(ranking_vip)
  
}