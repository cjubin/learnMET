#' Predict method for xgb_reg objects.
#'
#' Obtains predictions from a training/test split object, pre-processed
#' with [xgb_reg()].
#' 
#' @description 
#' Fit a gradient boosted trees model on an object of class `xgb_reg`.
#' Three hyperparameters (number of iterations = number of trees ; tree depth ;
#' learning rate) are tuned using the training set via Bayesian 
#' optimization with 5-folds cross-validation (k-folds CV). A model is fitted on
#' the training set using the best hyperparameters and model performance is evaluated on the 
#' test set. 
#' 
#' @param object an object of class `xgb_reg`
#' 
#' @param seed \code{integer} Seed value. 
#' 
#' @param inner_cv_reps \code{integer} Number of repeats of the k-folds CV
#'   for hyperparameter optimization.
#'   
#' @param inner_cv_folds \code{integer} Number k in the k-folds CV used for
#'   hyperparameter optimization.
#'   
#' @return 
#' 
#' 
#' @rdname fit_cv_split
#' @export
fit_cv_split.xgb_reg <- function(object,
                                 seed,
                                 inner_cv_reps = 1,
                                 inner_cv_folds = 5,
                                 ...) {
  if (class(object) != "xgb_reg") {
    stop("The object must be an object of the class 'xgb_reg'")
  }
  
  training = object[['training']]
  test = object[['test']]
  rec = object[['rec']]
  trait = as.character(rec$var_info[rec$var_info$role == 'outcome', 'variable'])
  
  
  
  # Define the prediction model to use
  
  xgboost_model <-
    parsnip::boost_tree(
      mode = "regression",
      trees = tune(),
      tree_depth = tune(),
      learn_rate = tune()
    ) %>%
    set_engine("xgboost", objective = "reg:squarederror") %>%
    translate()
  
  # Three hyperparameters are tuned for XGBoost.
  
  grid_hyperparameters <- parameters(trees(),
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
  
  folds <-
    vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
  
  
  cat('Optimization of hyperparameters for one training set has started.\n')
  set.seed(seed)
  opt_res <- wf %>%
    tune_bayes(
      resamples = folds,
      param_info = grid_hyperparameters,
      iter = 15,
      initial = 10,
      #iter = 20,
      #initial = 10,
      metrics = yardstick::metric_set(rmse),
      control = tune::control_bayes(verbose = FALSE, no_improve = 10)
    )
  
  cat('Optimizing hyperparameters for this training set: done!\n')
  
  # Retain the best hyperparameters and update the workflow with these
  # hyperparameters
  
  best_params <- opt_res %>%
    tune::select_best("rmse")
  
  
  model_final <- finalize_workflow(wf,
                                   best_params)
  
  # Fit the model on the train dataset and predict the test dataset
  
  fitted_model <-
    model_final %>%
    fit(data = training)
  
  cat('Fitting the training set: done!\n')
  
  # Obtain the variable importance with the gain metric
  
  predictors <- model_final %>%
    fit(data = training) %>%
    pull_workflow_fit()
  
  predictors <- predictors$fit$feature_names
  
  
  variable_importance_vip <- model_final %>%
    fit(data = training) %>%
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
  
  predictions_test <-
    as.data.frame(fitted_model %>% predict(new_data = test) %>% bind_cols(test))
  
  cor_pred_obs <-
    fitted_model %>% predict(new_data = test) %>% bind_cols(test) %>%
    group_by(IDenv) %>% summarize(COR = cor(.pred, get(trait), method = 'pearson'))
  
  rmse_pred_obs <-
    fitted_model %>% predict(new_data = test) %>% bind_cols(test) %>%
    group_by(IDenv) %>% summarize(RMSE = sqrt(mean((get(
      trait
    ) - .pred) ^ 2)))
  
  # Apply the trained data recipe
  rec <- prep(rec,strings_as_factors = FALSE)
  train = bake(rec, training)
  test = bake(rec, test)
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'best_hyperparameters' = as.data.frame(best_params),
      'training_transformed' = as.data.frame(train),
      'test_transformed' = as.data.frame(test),
      'ranking_vip' = ranking_vip
    ),
    class = 'res_fitted_split'
  )
  
  
  
  
  return(res_fitted_split)
  
  
  
}
