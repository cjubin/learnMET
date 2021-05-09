#' Fitting a model on the training set and predicting the test set.
#'
#' @description
#'
#' Hyperparameter optimization
#'
#' @param
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'



fitting_train_test_split_kernel <-
  function(split,
           seed,
           inner_cv_reps = 2,
           inner_cv_folds = 5,
           ...) {
    training = split[['training']]
    test = split[['test']]
    
    rec_G = split[['rec_G']]
    rec_E = split[['rec_E']]
    rec_GE = split[['rec_GE']]
    
    # Some settings common for all kernels to be trained
    
    metric <- yardstick::metric_set(rmse)
    
    ctrl_res <- stacks::control_stack_resamples()
    
    # Inner CV
    
    set.seed(seed)
    folds <-
      rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    # Define the prediction model to use
    
    svm_spec <-
      parsnip::svm_rbf(cost = tune("cost"),
                       rbf_sigma = tune("sigma")) %>%
      parsnip::set_engine("kernlab") %>%
      parsnip::set_mode("regression")
    
    
    # Add recipe and model definition to a workflow, for each of the kernel
    
    svm_wflow_G <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec) %>%
      workflows::add_recipe(rec_G)
    
    svm_wflow_E <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec) %>%
      workflows::add_recipe(rec_E)
    
    svm_wflow_GE <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec) %>%
      workflows::add_recipe(rec_GE)
    
    # tune cost and sigma with the inner CV for each of the kernels
    
    set.seed(seed)
    svm_res_E <-
      tune::tune_grid(
        svm_wflow_E,
        resamples = folds,
        grid = 10,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                            save_workflow = FALSE)
      )
    
    set.seed(seed)
    svm_res_GE <-
      tune::tune_grid(
        svm_wflow_GE,
        resamples = folds,
        grid = 10,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                            save_workflow = FALSE)
      )
    
    set.seed(seed)
    svm_res_G <-
      tune::tune_grid(
        svm_wflow_G,
        resamples = folds,
        grid = 10,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                            save_workflow = FALSE)
      )
    
    # Initialize a data stack using the stacks() function.
    
    METData_data_st <-
      stacks::stacks() %>%
      stacks::add_candidates(svm_res_G) %>%
      stacks::add_candidates(svm_res_E) %>%
      stacks::add_candidates(svm_res_GE)
    
    # Fit the stack
    
    METData_model_st <-
      METData_data_st %>%
      stacks::blend_predictions()
    
    METData_model_st <-
      METData_model_st %>%
      stacks::fit_members()
    
    
    predictions_test <-
      as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      sqrt(mean((predictions_test[, trait] - predictions_test[, '.pred']) ^ 2))
    
    
    return
    
    
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