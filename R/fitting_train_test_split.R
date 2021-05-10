#' Fitting a model on the training set and predicting the test set.
#'
#' @description
#'
#' Function creating a workflow based on the prediction model chosen (all 
#' methods other than multi-kernel learning), with hyperparameter optimization 
#' using a Bayesian approach on the training set. Once the best hyperparameters 
#' are identified via resampling, the model is fitted on the complete dataset. 
#' Finally the test set is predicted.
#'
#' @param split. A \code{split_processed} object containing:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set}
#'     \item{test}{\code{data.frame} Test set}
#'     \item{rec}{\code{recipe} object with the different steps to implement
#'      as processing steps on the training dataset. Same transformations applied
#'      on the test set.}
#'   }
#'     
#' @param prediction_method \code{character} Prediction method other than 
#'   kernel-based methods to use on the predictors defined by the recipe item of
#'   the split object.Options are 'xgboost' for XGBoost, 'rf' for Random Forest.
#'   Default is `xgboost`.  
#'   
#' @param seed \code{integer} Seed value.
#' 
#' @param inner_cv_reps \code{integer} Number of times to repeat the k-fold 
#'   partitioning used for the inner cross-validation for estimation of the best 
#'   hyperparameters. Default is 2.
#' 
#' @param inner_cv_folds \code{integer} Number of partitions of the training set
#'   for the inner cross-validation for estimation of the best hyperparameters.
#'   Default is 5.
#'   
#' @return a \code{list} object containing:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set.}
#'     \item{test}{\code{data.frame} Test set.}
#'     \item{predictions_df}{\code{data.frame} with original test dataset with 
#'      extra column containing predicted values.}
#'     \item{cor_pred_obs}{\code{numeric} Pearson's correlation between predicted
#'      and observed values of the test set.}
#'     \item{rmse_pred_obs}{\code{numeric} root mean square error between predicted and
#'      observed values of the test set.}
#'     \item{best_hyperparameters}{\code{tbl_df} the tuning parameter combination with
#'      the best performance values which was used to fit the final model on the
#'      training set.}
#'   }
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'



fitting_train_test_split <-
  function(split,
           prediction_method = c('xgboost'),
           seed,
           inner_cv_reps = 2,
           inner_cv_folds = 5,
           ...) {
    
    training = split[[1]]
    test = split[[2]]
    rec = split[[3]]
    
    
    
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
    }
    
    
  
    
    # Define folds for inner CV for optimization of hyperparameters (only used
    # on the training set)
    
    set.seed(seed)
    
    folds <-
      vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    
    cat('Optimization of hyperparameters for one training set has started.\n')
    set.seed(params$seed)
    opt_res <- wf %>%
      tune_bayes(
        resamples = folds,
        param_info = grid_hyperparameters,
        iter = 20,
        initial = 10,
        metrics = yardstick::metric_set(rmse),
        control = tune::control_bayes(verbose = FALSE, no_improve = 14)
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
    
    predictions_test <-
      as.data.frame(fitted_model %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      sqrt(mean((predictions_test[, trait] - predictions_test[, '.pred']) ^ 2))
    
    
    return(
      list(
        'training' = training,
        'test' = test,
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'best_hyperparameters' = best_params
        
      )
    )
    
    
    
  }