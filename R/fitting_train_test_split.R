#' Fitting a model on the training set and predicting the test set.
#' 
#' @details 
#' 
#' Hyperparameter optimization
#'
#' @param
#'
#'
#'
#'



fitting_train_test_split <- function(split, prediction_method,seed) {
  
  training = split[[1]]
  test = split[[2]]
  rec = split[[3]]
  
  
  prep(rec, training)
  
  
  if (prediction_method == 'xgboost') {
    # Define the prediction model to use
    
    xgboost_model <-
      parsnip::boost_tree(
        mode = "regression",
        trees = tune(),
        tree_depth = tune(),
        learn_rate = tune()
      ) %>%
      set_engine("xgboost", objective = "reg:linear") %>%
      translate()
    
    # Three hyperparameters are tuned for XGBoost.
    
    xgb_grid <- parameters(trees(),
                           learn_rate(),
                           tree_depth()) %>% update(
                             trees = trees(c(500, 4000)),
                             learn_rate = learn_rate(range(c(5e-4, 0.05)), trans = NULL),
                             tree_depth = tree_depth(c(2, 20))
                           )
    
    # Workflow with recipe
    
    wf <- workflow() %>%
      add_model(xgboost_model) %>%
      add_recipe(rec)
    
    # Define folds for inner CV for optimization of hyperparameters (only used
    # on the training set)
    set.seed(seed)
    
    folds <- vfold_cv(training, repeats = 1, v = 5)
    
    # 
    cat('Optimization of hyperparameters for one training set has started.\n')
    
    opt_res <- wf %>%
      tune_bayes(
        resamples = folds,
        param_info = xgb_grid,
        iter = 5,
        initial = 8,
        metrics = yardstick::metric_set(rmse),
        control = tune::control_bayes(verbose = FALSE, no_improve = 5)
      )
    
    cat('Optimizing hyperparameters for this training set: done!\n')
    
    # Retain the best hyperparameters and update the workflow with these
    # hyperparameters
    
    xgboost_best_params <- opt_res %>%
      tune::select_best("rmse")
    
    
    xgboost_model_final <- finalize_workflow(wf,
                                             xgboost_best_params)
    
    # Fit the model on the train dataset and predict the test dataset
    
    fitted_model <-
      xgboost_model_final %>%
      fit(data = training)
    
    cat('Fitting the training set: done!\n')
    
    predictions_test <-
      as.data.frame(fitted_model %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      sqrt(mean((
        predictions_test[, trait] - predictions_test[, '.pred']
      ) ^ 2))
    
    
    return(
      list(
        'training' = training,
        'test' = test,
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs
      )
    )
  }
  
  
}