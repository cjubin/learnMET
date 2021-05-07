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
    
    training_G = split[['training_G']]
    training_E = split[['training_E']]
    training_GE = split[['training_GE']]
    test_G = split[['test_G']]
    test_E = split[['test_E']]
    test_GE = split[['test_GE']]
    rec_G = split[['rec_G']]
    rec_E = split[['rec_E']]
    rec_GE = split[['rec_GE']]
    
    # Some settings common for all kernels to be trained
    
    metric <- yardstick::metric_set(rmse)
    ctrl_grid <- stacks::control_stack_grid()
    ctrl_res <- stacks::control_stack_resamples()
    
    # Inner CV
    
    set.seed(seed)
    folds_G <- rsample::vfold_cv(training_G, repeats = inner_cv_reps, v = inner_cv_folds)
    set.seed(seed)
    folds_E <- rsample::vfold_cv(training_E, repeats = inner_cv_reps, v = inner_cv_folds)
    set.seed(seed)
    folds_GE <- rsample::vfold_cv(training_GE, repeats = inner_cv_reps, v = inner_cv_folds)
    
    
    # Define the prediction model to use
    
    svm_spec <- 
      parsnip::svm_rbf(
        cost = tune("cost"), 
        rbf_sigma = tune("sigma")
      ) %>%
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
    svm_res_G <- 
      tune::tune_grid(
        svm_wflow_G, 
        resamples = folds_G, 
        grid = 6,
        metrics = metric,
        control = ctrl_grid
      )
      
    set.seed(seed)
    svm_res_E <- 
      tune::tune_grid(
        svm_wflow_E, 
        resamples = folds_E, 
        grid = 6,
        metrics = metric,
        control = ctrl_grid
      )
    
    set.seed(seed)
    svm_res_GE <- 
      tune::tune_grid(
        svm_wflow_GE, 
        resamples = folds_GE, 
        grid = 6,
        metrics = metric,
        control = ctrl_grid
      )
    
    
    # Initialize a data stack using the stacks() function.
    
    METData_data_st <- 
      stacks() %>%
      add_candidates(svm_res_G) %>%
      add_candidates(svm_res_E) %>%
      add_candidates(svm_res_GE)
    
    # Fit the stack
    
    METData_model_st <-
      METData_data_st %>%
      blend_predictions()
    
    METData_model_st <-
      METData_model_st %>%
      fit_members()
    
    
    tree_frogs_test <- 
      tree_frogs_test %>%
      bind_cols(predict(tree_frogs_model_st, .))
    
    # Define folds for inner CV for optimization of hyperparameters (only used
    # on the training set)
    
    set.seed(seed)
    
    folds <-
      vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    #
    cat('Optimization of hyperparameters for one training set has started.\n')
    
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