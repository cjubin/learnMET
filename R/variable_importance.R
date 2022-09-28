#' Compute variable importance according to the machine learning algorithm used
#'
#' @description
#' Variable importance can be calculated based on model-specific and
#' model-agnostic approaches
#'
#'
#' @name variable_importance_split
#'
#' @param object an object of class `res_fitted_split`
#'
#' @param type `model_specific` or `model_agnostic`
#'
#' @param permutations By default, equal to 10.
#'
#' @param unseen_data
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @references
#'
#' \insertRef{breiman2001random}{learnMET}
#' \insertRef{molnar2022}{learnMET}
#' @export
variable_importance_split <- function(object, ...) {
  UseMethod("variable_importance_split")
}


#' @rdname variable_importance_split
#' @export
variable_importance_split.default <- function(object, ...) {
  stop('not implemented')
  
}

#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_DL_reg_1 <-
  function(object) {
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
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
      label = "DL_model_vip",
      verbose = FALSE
    )
    
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        loss_function = DALEX::loss_root_mean_square,
        type = 'difference',
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
variable_importance_split.fitted_DL_reg_2 <-
  function(object) {
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
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
      label = "DL_model_vip",
      verbose = FALSE
    )
    
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        loss_function = DALEX::loss_root_mean_square,
        type = 'difference',
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
variable_importance_split.fitted_xgb_reg_1 <-
  function(object,
           path_plot,
           type = 'model_specific',
           permutations = 10,
           unseen_data = F) {
    if (type == 'model_specific') {
      # Obtain the variable importance with the gain metric
      cat('Variable importance with gain metric\n')
      
      model <- fitted_split$fitted_model
      trait <-
        as.character(fitted_split$fitted_model$pre$actions$recipe$recipe$var_info[fitted_split$fitted_model$pre$actions$recipe$recipe$var_info$role ==
                                                                                    'outcome', 'variable'])
      y_train <- as.matrix(fitted_split$training %>%
                             dplyr::select(all_of(trait)))
      x_train <- fitted_split$training
      
      
      predictors <- model %>%
        parsnip::fit(data = x_train) %>%
        workflows::pull_workflow_fit()
      
      predictors <- predictors$fit$feature_names
      
      
      ranking_vip <- as.data.frame(
        model %>%
          parsnip::fit(data = x_train) %>%
          workflows::pull_workflow_fit() %>%
          vip::vi()
      )
      
      
      if (length(predictors[which(predictors %notin% ranking_vip$Variable)])
          >
          0) {
        remaining <-
          cbind(as.vector(predictors[which(predictors %notin% ranking_vip$Variable)]), as.numeric(0))
        colnames(remaining) <- colnames(ranking_vip)
        ranking_vip <- rbind(ranking_vip, remaining)
      }
      ranking_vip$Importance <- as.numeric(ranking_vip$Importance)
      
      colnames(ranking_vip) <- c('Variable', 'Importance')
      
      plot_results_vip(x = ranking_vip,
                       path_plot = path_plot,
                       type = type)
      
      
      return(ranking_vip)
      
    }
    
    if (type == 'model_agnostic') {
      cat('Permutation feature importance - Nb of permutations: ',permutations,'\n')
      
      model <- fitted_split$fitted_model
      trait <-
        colnames(workflows::extract_mold(fitted_split$fitted_model)$outcome)
      
      predictors <-
        colnames(workflows::extract_mold(fitted_split$fitted_model)$predictor)
      
      if (unseen_data) {
        # use test set if permutation feature importance evaluated on test set
        x <- fitted_split$test
        
        y <- as.numeric(as.data.frame(fitted_split$test %>%
                                        dplyr::select(all_of(trait)))[, 1])
        
      }  else{
        # otherwise use the training set with which the model was fitted to
        y <- as.numeric(as.data.frame(fitted_split$training %>%
                                        dplyr::select(all_of(trait)))[, 1])
        x <- fitted_split$training
      }
      
      # Permutation VIP function
      res_permutations <- permutation_based_vip(
        model,
        x = x,
        y = y,
        permutations = permutations,
        predictors = predictors,
        path_plot = path_plot
      )
      
      
      return(res_permutations)
    }
  }

#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_xgb_reg_2 <-
  function(object) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
    
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
  function(object) {
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- as.numeric(fitted_split$y_train[, trait])
    x_train <- fitted_split$x_train
    env_predictors <- fitted_split$env_predictors
    
    print(
      'Variable importance (permutation-based) will only be computed for environmental features.'
    )
    vars = env_predictors
    
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(as.numeric(as.vector(as.data.frame(results)[, '.pred'])))
    }
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = y_train,
      predict_function = pred_wrapper,
      label = "stacked_model_vip",
      verbose = FALSE
    )
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        variables = vars,
        loss_function = DALEX::loss_root_mean_square,
        type = 'difference',
        B = 10
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

#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_stacking_reg_2 <-
  function(object) {
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- as.numeric(fitted_split$y_train[, trait])
    x_train <- fitted_split$x_train
    env_predictors <- fitted_split$env_predictors
    
    print(
      'Variable importance (permutation-based) will only be computed for environmental features.'
    )
    vars = env_predictors
    
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(as.numeric(as.vector(as.data.frame(results)[, '.pred'])))
    }
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = y_train,
      predict_function = pred_wrapper,
      label = "stacked_model_vip",
      verbose = FALSE
    )
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        variables = vars,
        loss_function = DALEX::loss_root_mean_square,
        type = 'difference',
        B = 10
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


#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_stacking_reg_3 <-
  function(object) {
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- as.numeric(fitted_split$y_train[, trait])
    x_train <- fitted_split$x_train
    env_predictors <- fitted_split$env_predictors
    
    print(
      'Variable importance (permutation-based) will only be computed for environmental features.'
    )
    vars = env_predictors
    
    # create custom predict function
    pred_wrapper <- function(model, newdata)  {
      results <- model %>% predict(new_data = newdata) %>%
        as.vector()
      
      return(as.numeric(as.vector(as.data.frame(results)[, '.pred'])))
    }
    
    explainer <- DALEX::explain(
      model = model,
      data = x_train,
      y = y_train,
      predict_function = pred_wrapper,
      label = "stacked_model_vip",
      verbose = FALSE
    )
    ranking_vip <-
      DALEX::model_parts(
        explainer,
        N = NULL,
        variables = vars,
        loss_function = DALEX::loss_root_mean_square,
        B = 10,
        type = 'difference'
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

#' @rdname variable_importance_split
#' @export
variable_importance_split.fitted_rf_reg_1 <-
  function(object) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
    
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
variable_importance_split.fitted_rf_reg_2 <-
  function(object) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
    
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
variable_importance_split.fitted_rf_reg_2 <-
  function(object) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
    
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
variable_importance_split.fitted_rf_reg_1 <-
  function(object) {
    # Obtain the variable importance with the gain metric
    
    model <- fitted_split$fitted_model
    trait <- fitted_split$trait
    y_train <- fitted_split$y_train
    x_train <- fitted_split$x_train
    
    
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
