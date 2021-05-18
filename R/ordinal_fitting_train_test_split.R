#' Fitting a model on the training set and predicting the test set.
#'
#' @description
#' Function creating a workflow based on the prediction model chosen (all
#' methods other than multi-kernel learning), with hyperparameter optimization
#' using a Bayesian approach on the training set. Once the best hyperparameters
#' are identified via resampling, the model is fitted on the complete dataset.
#' Finally the test set is predicted.
#'
#' @param split An object of class \code{split_processed} object with the
#'   following items:
#'   * **training**: \code{data.frame} Training set.
#'   * **test**: \code{data.frame} Test set.
#'   * **rec**: \code{recipe} object with the different steps to implement
#'   as processing steps on the training dataset. Same transformations applied
#'   on the test set.
#'   }
#'
#' @param prediction_method \code{character} Prediction method other than
#'   kernel-based methods to use on the predictors defined by the recipe item of
#'   the split object.Options are 'xgboost' for XGBoost, 'rf' for Random Forest.
#'   Default is `xgboost`.
#'
#' @param treatment_ordinal_data 'naive_multiclass','ordinal_prediction',
#' 'regression'
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
#' @param return_finalized_train_test_sets a \code{logical} whether the trained
#'   dataset and the test dataset resulting from preprocessing operations
#'   should be returned in the final `res_fitted_split` object. Default is
#'   `FALSE`.
#'
#' @return a \code{list} object of class \code{res_fitted_split} with the
#'   following items:
#'   * **predictions_df**: \code{data.frame} with original test dataset with
#'      extra column containing predicted values.
#'   * **cor_pred_obs**: \code{numeric} Pearson's correlation between predicted
#'      and observed values of the test set.
#'   * **rmse_pred_obs**: \code{numeric} root mean square error between
#'   predicted and observed values of the test set.
#'   * **best_hyperparameters**: a \code{tbl_df} giving the tuning parameter
#'   combination with the best performance values which was used to fit the
#'   final model on the training set.
#'
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'



ordinal_fitting_train_test_split <-
  function(split,
           prediction_method = c('xgboost'),
           treatment_ordinal_data = c('ordinal_prediction'),
           seed,
           nb_classes,
           inner_cv_reps = 2,
           inner_cv_folds = 5,
           return_finalized_train_test_sets = F,
           ...) {
    training = split[[1]]
    test = split[[2]]
    rec = split[[3]]
    
    
    if (prediction_method == 'xgboost') {
      # Define the prediction model to use
      if (treatment_ordinal_data == 'ordinal_prediction') {
        xgboost_model <-
          parsnip::boost_tree(
            mode = "classification",
            trees = tune(),
            tree_depth = tune(),
            learn_rate = tune()
          ) %>%
          set_engine("xgboost", objective = "binary:logistic") %>%
          translate()
      } else if (treatment_ordinal_data == 'naive_multiclass') {
        xgboost_model <-
          parsnip::boost_tree(
            mode = "classification",
            trees = tune(),
            tree_depth = tune(),
            learn_rate = tune()
          ) %>%
          set_engine(
            "xgboost",
            objective = "multi:softprob",
            lambda = 0,
            alpha = 1,
            num_class = nb_classes
          ) %>%
          translate()
      } else {
        xgboost_model <-
          parsnip::boost_tree(
            mode = "regression",
            trees = tune(),
            tree_depth = tune(),
            learn_rate = tune()
          ) %>%
          set_engine(
            "xgboost",
            objective = "reg:squarederror"
          ) %>%
          translate()
        
      }
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
        add_formula(trait_binary ~ .)
    }
    
    
    if (treatment_ordinal_data == 'ordinal_prediction') {
      # For each k - 1 ordinal value we will fit a binary classification problem
      
      binary_class_k <- function(k, training, test, wf) {
        rec <- rec %>%
          step_mutate(trait_binary = case_when(get(trait) <= k ~ 0,
                                               get(trait) > k ~ 1),
                      role = 'outcome') %>%
          step_rm(all_of(trait)) %>%
          step_rm(geno_ID) %>%
          step_rm(IDenv) %>%
          step_mutate(trait_binary = as.factor(trait_binary))
        
        prepped <- prep(rec)
        training <- bake(prepped, training)
        test <- bake(prepped, test)
        
        
        
        
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
            iter = 2,
            initial = 4,
            #iter = 20,
            #initial = 10,
            metrics = yardstick::metric_set(
              recall,
              precision,
              f_meas,
              accuracy,
              kap,
              roc_auc,
              sens,
              spec
            ),
            control = tune::control_bayes(verbose = FALSE, no_improve = 14)
          )
        
        cat('Optimizing hyperparameters for this training set: done!\n')
        
        # Retain the best hyperparameters and update the workflow with these
        # hyperparameters
        
        best_params <- opt_res %>%
          tune::select_best("roc_auc")
        
        
        model_final <- finalize_workflow(wf,
                                         best_params)
        
        # Fit the model on the train dataset and predict the test dataset
        
        fitted_model <-
          model_final %>%
          fit(data = training)
        
        cat('Fitting the training set: done!\n')
        
        predictions_test <-
          as.data.frame(fitted_model %>% predict(new_data = test, type = "prob") %>% bind_cols(test))
        
        colnames(predictions_test)[which(colnames(predictions_test) == ".pred_0")] <-
          paste0('.pred_', 'inf_or_equal_to_', k)
        colnames(predictions_test)[which(colnames(predictions_test) == ".pred_1")] <-
          paste0('.pred_', 'sup_to_', k)
        
        return(predictions_test[, c(1, 2)])
        
      }
      
      res_all_binary_classifiers <-
        lapply(seq_len(nb_classes-1), function(x) {
          binary_class_k(
            k = x,
            training = training,
            test = test,
            wf = wf
          )
        })
      
      
      
      
      
      
      # Return final list of class res_fitted_split
      if (return_finalized_train_test_sets) {
        res_fitted_split <- list(
          'parameters_collection_G' = parameters_collection_G,
          'parameters_collection_E' = parameters_collection_E,
          'parameters_collection_GE' = parameters_collection_GE,
          'predictions_df' = predictions_test,
          'cor_pred_obs' = cor_pred_obs,
          'rmse_pred_obs' = rmse_pred_obs,
          'training' = train,
          'test' = test
          
          
        )
      } else{
        res_fitted_split <- list(
          'parameters_collection_G' = parameters_collection_G,
          'parameters_collection_E' = parameters_collection_E,
          'parameters_collection_GE' = parameters_collection_GE,
          'predictions_df' = predictions_test,
          'cor_pred_obs' = cor_pred_obs,
          'rmse_pred_obs' = rmse_pred_obs
          
        )
        
      }
      class(res_fitted_split) <- c('res_fitted_split', 'list')
      return(res_fitted_split)
      
      
      
      
      
    }
  }
