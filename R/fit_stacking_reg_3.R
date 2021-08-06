#' @rdname fit_cv_split
#' @export
fit_cv_split.stacking_reg_3 <- function (object,
                                           seed,
                                           inner_cv_reps = 1,
                                           inner_cv_folds = 5,
                                           kernel_G = 'rbf',
                                           kernel_E = 'rbf',
                                           ...) {
  # Case if the GE kernel was built (length = 5)
 
  
    training = object[['training']]
    test = object[['test']]
    
    rec_G = object[['rec_G']]
    rec_E = object[['rec_E']]
    rec_all = object[['rec_all']]
    
    trait = as.character(rec_G$var_info[rec_G$var_info$role == 'outcome', 'variable'])
    
    # Some settings common for all kernels to be trained
    
    metric <- yardstick::metric_set(rmse)
    
    ctrl_res <- stacks::control_stack_resamples()
    
    # Inner CV
    
    set.seed(seed)
    folds <-
      rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    # Define the prediction model to use: support vector machine with SNPs,
    # support vector machine with environmental covariables, and PLS with
    # SNPs and ECs.
    # Define the space-filling design for grid search.
    
    pls_all_predictors <-
      plsmod::pls(num_comp = tune()) %>% 
      parsnip::set_engine("mixOmics") %>% 
      parsnip::set_mode("regression")
    
    grid_model_pls <- tibble(num_comp = seq(from = 1, to = 20, by = 2))

      
    svm_spec_rbf <-
      parsnip::svm_rbf(cost = tune("cost"),
                       rbf_sigma = tune("sigma")) %>%
      parsnip::set_engine("kernlab") %>%
      parsnip::set_mode("regression")
    
    svm_spec_polynomial <-
      parsnip::svm_poly(
        cost = tune("cost"),
        degree = tune("degree"),
        scale_factor = tune('scale_factor')
      ) %>%
      parsnip::set_engine("kernlab") %>%
      parsnip::set_mode("regression")
    
    svm_spec_linear <-
      parsnip::svm_linear(cost = tune("cost")) %>%
      parsnip::set_engine("LiblineaR") %>%
      parsnip::set_mode("regression")
    
    
    if (kernel_G == 'rbf') {
      svm_spec_G <- svm_spec_rbf
      grid_model_G <- 8
    } else if (kernel_G == 'polynomial') {
      svm_spec_G <- svm_spec_polynomial
      grid_model_G <- 14
    } else{
      svm_spec_G <- svm_spec_linear
      grid_model_G <- 6
    }
    
    if (kernel_E == 'rbf') {
      svm_spec_E <- svm_spec_rbf
      grid_model_E <- 8
    } else if (kernel_E == 'polynomial') {
      svm_spec_E <- svm_spec_polynomial
      grid_model_E <- 14
    } else{
      svm_spec_E <- svm_spec_linear
      grid_model_E <- 6
    }
    
    
    
    # Add recipe and model definition to a workflow, for each of the kernel
    
    pls_wflow_all <-
      workflows::workflow() %>%
      workflows::add_model(pls_all_predictors) %>%
      workflows::add_recipe(rec_all)
    
    svm_wflow_G <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_G) %>%
      workflows::add_recipe(rec_G)
    
    svm_wflow_E <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_E) %>%
      workflows::add_recipe(rec_E)
    
    
    # tune cost and sigma with the inner CV for each of the kernels
    print(grid_model_pls)
    set.seed(seed)
    pls_res_all <-
      tune::tune_grid(
        pls_wflow_all,
        resamples = folds,
        grid = grid_model_pls,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    cat('PLS with all predictors done!')
    
    
    set.seed(seed)
    svm_res_E <-
      tune::tune_grid(
        svm_wflow_E,
        resamples = folds,
        grid = grid_model_G,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    cat('Support vector regression with env. kernel done!')
    
    
    set.seed(seed)
    svm_res_G <-
      tune::tune_grid(
        svm_wflow_G,
        resamples = folds,
        grid = grid_model_G,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    cat('Support vector regression with genomic kernel done!')
    
    
    # Initialize a data stack using the stacks() function.
    
    METData_data_st <-
      stacks::stacks() %>%
      stacks::add_candidates(svm_res_G) %>%
      stacks::add_candidates(svm_res_E) %>%
      stacks::add_candidates(pls_res_all)
    
    # Fit the stack
    
    METData_model_st <-
      METData_data_st %>%
      stacks::blend_predictions()
    
    METData_model_st <-
      METData_model_st %>%
      stacks::fit_members()
    
    # Identify which model configurations were assigned what stacking coefficients
    parameters_collection_G <-
      METData_model_st %>% stacks::collect_parameters('svm_res_G')
    parameters_collection_E <-
      METData_model_st %>% stacks::collect_parameters('svm_res_E')
    parameters_collection_all_predictors <-
      METData_model_st %>% stacks::collect_parameters('pls_res_all')
    
    
    # Predictions and metrics calculated on a per-environment basis
    
    predictions_test <-
      as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      METData_model_st %>% predict(new_data = test) %>% bind_cols(test) %>%
      group_by(IDenv) %>% summarize(COR = cor(.pred, get(trait), method = 'pearson'))
    print(cor_pred_obs)
    
    #cor_pred_obs <-
    #  cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      METData_model_st %>% predict(new_data = test) %>% bind_cols(test) %>%
      group_by(IDenv) %>% summarize(RMSE = sqrt(mean((
        get(trait) - .pred
      ) ^ 2)))
    
    # Apply the trained data recipe
    rec_all <- prep(rec_all)
    train_all = bake(rec_all, training)
    test_all = bake(rec_all, test)
    
    
    # Return final list of class res_fitted_split
    
    res_fitted_split <- structure(
      list(
        'prediction_method' = class(object),
        'parameters_collection_G' = as.data.frame(parameters_collection_G),
        'parameters_collection_E' = as.data.frame(parameters_collection_E),
        'parameters_collection_all_predictors' = as.data.frame(parameters_collection_all_predictors),
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'training_all_transformed' = as.data.frame(train_all),
        'test_all_transformed' = as.data.frame(test_all)
      ),
      class = 'res_fitted_split'
    )
  
  
  return(res_fitted_split)
  
  
}
