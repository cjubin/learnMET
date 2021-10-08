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
variable_importance_split.fitted_DL_reg <-
  function(fitted_obj_for_vip) {
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
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = as.numeric(y_train),
      predict_function = pred_wrapper,
      type = 'regression',
      label = "DL_model_vip"
    )
    
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        loss_function = DALEX::loss_root_mean_square,
        B = 10
      )
    ranking_vip_avg <- as.data.frame(
      as.data.frame(ranking_vip) %>%
        group_by(variable) %>%
        dplyr::summarize(Importance = mean(dropout_loss, na.rm = TRUE))
    )
    colnames(ranking_vip_avg) <- c('Variable', 'Importance')
    return(ranking_vip_avg)
  }

#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_xgb_reg <-
  function(fitted_obj_for_vip) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_obj_for_vip$model
    trait <- fitted_obj_for_vip$trait
    y_train <- fitted_obj_for_vip$y_train
    x_train <- fitted_obj_for_vip$x_train
    
    
    predictors <- model %>%
      fit(data = x_train) %>%
      pull_workflow_fit()
    
    predictors <- predictors$fit$feature_names
    
    
    ranking_vip <- as.data.frame(model %>%
                                   fit(data = x_train) %>%
                                   pull_workflow_fit() %>%
                                   vip::vi(method = 'model'))
    
    
    
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
variable_importance_split.fitted_stacking_reg_1 <-
  function(fitted_obj_for_vip) {
    model <- fitted_obj_for_vip$model
    trait <- fitted_obj_for_vip$trait
    y_train <- as.numeric(fitted_obj_for_vip$y_train[, trait])
    x_train <- fitted_obj_for_vip$x_train
    env_predictors <- fitted_obj_for_vip$env_predictors
    
    print('Variable importance (permutation-based) will only be computed for environmental features.')
    vars = env_predictors
      
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(results[, '.pred'])
    }
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = y_train,
      predict_function = pred_wrapper,
      label = "stacked_model_vip"
    )
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        variables=vars,
        loss_function = DALEX::loss_root_mean_square,
        B = 1
      )
    ranking_vip_avg <- as.data.frame(
      as.data.frame(ranking_vip) %>%
        group_by(variable) %>%
        dplyr::summarize(Importance = mean(dropout_loss, na.rm =
                                             TRUE))
    )
    colnames(ranking_vip_avg) <- c('Variable', 'Importance')
    return(ranking_vip_avg)
    
    return(ranking_vip)
    
  }