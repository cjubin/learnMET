#' Processing and selecting predictors to fit the model
#'
#' @param
#'
#'
#'
#'



fitting_train_test_split <- function(split, prediction_method) {
  
  training = split[[1]]
  test = split[[2]]
  rec = split[[3]]
  
  
  prep(rec,training)
  
  
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
                             trees = trees(c(300, 1500)),
                             learn_rate = learn_rate(range(c(0.005, 0.05)), trans = NULL),
                             tree_depth = tree_depth(c(2, 20))
                           )
    
    # Workflow with recipe
    
    wf <- workflow() %>%
      add_model(xgboost_model) %>%
      add_recipe(rec)
    
    # Define folds for inner CV for optimization of hyperparameters (only used 
    # on the training set)
    
    folds <- vfold_cv(training, repeats = 1, v = 5)
    
    opt_res <- wf %>%
      tune_bayes(
        resamples = folds,
        param_info = xgb_grid,
        iter = 6,
        initial = 8,
        metrics = yardstick::metric_set(rmse),
        control = tune::control_bayes(verbose = TRUE, no_improve = 5)
      )
    
    # Retain the best hyperparameters and update the workflow with these
    # hyperparameters
    
    xgboost_best_params <- opt_res %>%
      tune::select_best("rmse")

    
    xgboost_model_final <- finalize_workflow(
      wf,
      xgboost_best_params
    )
    
    # Fit the model on the train dataset and predict the test dataset
    
    fitted_model <- 
      xgboost_model_final %>%
      fit(data = training) 
    
    predictions_test <- fitted_model %>% predict(new_data = test)
    
    cor_pred_obs <- cor(predictions_test$.pred,test[,trait],method = 'pearson')
    
    return(list(training,test,cor_pred_obs))
  }
  
  
}