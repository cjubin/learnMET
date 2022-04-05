#' @rdname fit_cv_split
#' @export
fit_cv_split.stacking_reg_3 <- function (object,
                                         seed,
                                         inner_cv_reps = 1,
                                         inner_cv_folds = 4,
                                         kernel_E = 'polynomial',
                                         save_model = F,
                                         penalty = 1, #lambda
                                         mixture = 0, #alpha
                                         ...) {
  
  cat('Kernel for E is', kernel_E,'\n')
  
  
  training = object[['training']]
  test = object[['test']]
  
  rec_G = object[['rec_G']]
  rec_E = object[['rec_E']]
  rec_GE = object[['rec_GE']]
  
  trait = as.character(rec_G$var_info[rec_G$var_info$role == 'outcome', 'variable'])
  
  env_predictors = colnames(
    recipes::bake(recipes::prep(rec_E), new_data = training) %>%
      dplyr::select(-IDenv, -tidyselect::all_of(trait))
  )
  
  
  
  # Some settings common for all kernels to be trained
  
  metric <- yardstick::metric_set(yardstick::rmse)
  
  ctrl_res <- stacks::control_stack_resamples()
  
  # Inner CV
  
  set.seed(seed)
  folds <-
    rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
  
  # Define the prediction model to use: support vector machine with 3 different
  # subset features.
  # Define the space-filling design for grid search.
  
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
  
  xgboost_model <-
    parsnip::boost_tree(
      mode = "regression",
      trees = 1000,
      tree_depth = tune(),
      learn_rate = tune()
    ) %>%
    parsnip::set_engine("xgboost", objective = "reg:squarederror") %>%
    parsnip::translate()
  
  en_model <- parsnip::linear_reg(penalty = tune(), mixture = tune()) %>%  #elastic net
    parsnip::set_engine("glmnet")
  
  mixture_param <- tune::parameters(dials::penalty(), dials::mixture())
  
  grid_model_G <- dials::grid_max_entropy(mixture_param, size = 10)
  
   
  
  
  if (kernel_E == 'rbf') {
    svm_spec_E <- svm_spec_rbf
    grid_model_E <- 5
  } else if (kernel_E == 'polynomial') {
    svm_spec_E <- svm_spec_polynomial
    grid_model_E <- 5
  } else{
    svm_spec_E <- svm_spec_linear
    grid_model_E <- 5
  }
  
  # grid specification for xgb
  xgboost_params <- 
    tune::parameters(dials::learn_rate(),
                     dials::tree_depth()) %>% update(
                       learn_rate = dials::learn_rate(range(c(5e-4, 0.05)), trans = NULL),
                       tree_depth = dials::tree_depth(c(2, 12))
                     )
  
  xgboost_grid <- 
    dials::grid_max_entropy(
      xgboost_params, 
      size = 8
    )
  
  
  
  # Add recipe and model definition to a workflow, for each of the kernel
  
  en_wflow_G <-
    workflows::workflow() %>%
    workflows::add_model(en_model) %>%
    workflows::add_recipe(rec_G)
  
  svm_wflow_E <-
    workflows::workflow() %>%
    workflows::add_model(svm_spec_E) %>%
    workflows::add_recipe(rec_E)
  
  xgb_wflow_GE <-
    workflows::workflow() %>%
    workflows::add_model(xgboost_model) %>%
    workflows::add_recipe(rec_GE)
  
  # tune cost and sigma with the inner CV for each of the kernels
  
  set.seed(seed)
  
  set.seed(seed)
  svm_res_E <-
    tune::tune_grid(
      svm_wflow_E,
      resamples = folds,
      grid = grid_model_E,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with env. kernel done!')
  
  
  set.seed(seed)
  en_res_G <-
    tune::tune_grid(
      en_wflow_G,
      resamples = folds,
      grid = grid_model_G,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with genomic kernel done!')
  
  
  
  set.seed(seed)
  xgb_res_GE <-
    tune::tune_grid(
      xgb_wflow_GE,
      resamples = folds,
      grid = xgboost_grid,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with GxE kernel done!')
  
  # Initialize a data stack using the stacks() function.
  
  METData_data_st <-
    stacks::stacks() %>%
    stacks::add_candidates(en_res_G) %>%
    stacks::add_candidates(svm_res_E) %>%
    stacks::add_candidates(xgb_res_GE)
  
  # Fit the stack
  
  METData_model_st <-
    METData_data_st %>%
    stacks::blend_predictions()
  
  METData_model_st <-
    METData_model_st %>%
    stacks::fit_members()
  
  # Identify which model configurations were assigned what stacking coefficients
  parameters_collection_G <-
    METData_model_st %>% stacks::collect_parameters('en_res_G')
  parameters_collection_E <-
    METData_model_st %>% stacks::collect_parameters('svm_res_E')
  parameters_collection_GE <-
    METData_model_st %>% stacks::collect_parameters('xgb_res_GE')
  
  
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
    METData_model_st %>% predict(new_data = test) %>% dplyr::bind_cols(test) %>%
    dplyr::group_by(IDenv) %>% dplyr::summarize(RMSE = sqrt(mean((get(
      trait
    ) - .pred) ^ 2)))
  
  
  # Apply the trained data recipe
  rec_GE <- recipes::prep(rec_GE)
  train_GE = recipes::bake(rec_GE, training)
  test_GE = recipes::bake(rec_GE, test)
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'training' = as.data.frame(training),
      'test' = as.data.frame(test)
    ),
    class = c('res_fitted_split','fitted_stacking_reg_3')
  )
  
  
  if (save_model) {
    res_fitted_split[["fitted_model"]] = METData_model_st
  }  else{
    res_fitted_split["fitted_model"] = list(NULL)
  }
  
  return(res_fitted_split)
  
  
  
  
  
  return(res_fitted_split)
  
  
}




#' @rdname fit_split
#' @export
fit_split.stacking_reg_3 <- function (object,
                                         seed,
                                         inner_cv_reps = 1,
                                         inner_cv_folds = 5,
                                         kernel_G = 'linear',
                                         kernel_E = 'polynomial',
                                         save_model = F,
                                         ...) {
 
  
  training = object[['training']]
  test = object[['test']]
  
  rec_G = object[['rec_G']]
  rec_E = object[['rec_E']]
  rec_GE = object[['rec_GE']]
  
  trait = as.character(rec_G$var_info[rec_G$var_info$role == 'outcome', 'variable'])
  
  env_predictors = colnames(
    recipes::bake(recipes::prep(rec_E), new_data = training) %>%
      dplyr::select(-IDenv, -tidyselect::all_of(trait))
  )
  
  
  
  # Some settings common for all kernels to be trained
  
  metric <- yardstick::metric_set(yardstick::rmse)
  
  ctrl_res <- stacks::control_stack_resamples()
  
  # Inner CV
  
  set.seed(seed)
  folds <-
    rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
  
  # Define the prediction model to use: support vector machine with 3 different
  # subset features.
  # Define the space-filling design for grid search.
  
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
  
  xgboost_model <-
    parsnip::boost_tree(
      mode = "regression",
      trees = 2000,
      tree_depth = tune(),
      learn_rate = tune()
    ) %>%
    parsnip::set_engine("xgboost", objective = "reg:squarederror") %>%
    parsnip::translate()
  
  
  
  en_model <- parsnip::linear_reg(penalty = tune(), mixture = tune()) %>%  #elastic net
    parsnip::set_engine("glmnet")
  
  mixture_param <- tune::parameters(dials::penalty(), dials::mixture())
  
  grid_model_G <- dials::grid_max_entropy(mixture_param, size = 10)
  
  
  
  
  if (kernel_E == 'rbf') {
    svm_spec_E <- svm_spec_rbf
    grid_model_E <- 5
  } else if (kernel_E == 'polynomial') {
    svm_spec_E <- svm_spec_polynomial
    grid_model_E <- 5
  } else{
    svm_spec_E <- svm_spec_linear
    grid_model_E <- 5
  }
  
  # grid specification for xgb
  xgboost_params <- 
    tune::parameters(dials::learn_rate(),
                     dials::tree_depth()) %>% update(
                       learn_rate = dials::learn_rate(range(c(5e-4, 0.05)), trans = NULL),
                       tree_depth = dials::tree_depth(c(2, 12))
                     )
  
  
  xgboost_grid <- 
    dials::grid_max_entropy(
      xgboost_params, 
      size = 6
    )
  
  
  
  # Add recipe and model definition to a workflow, for each of the kernel
  
  en_wflow_G <-
    workflows::workflow() %>%
    workflows::add_model(en_model) %>%
    workflows::add_recipe(rec_G)
  
  svm_wflow_E <-
    workflows::workflow() %>%
    workflows::add_model(svm_spec_E) %>%
    workflows::add_recipe(rec_E)
  
  xgb_wflow_GE <-
    workflows::workflow() %>%
    workflows::add_model(xgboost_model) %>%
    workflows::add_recipe(rec_GE)
  
  # tune cost and sigma with the inner CV for each of the kernels
  
  set.seed(seed)
  
  set.seed(seed)
  svm_res_E <-
    tune::tune_grid(
      svm_wflow_E,
      resamples = folds,
      grid = grid_model_E,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with env. kernel done!\n')
  
  
  set.seed(seed)
  en_res_G <-
    tune::tune_grid(
      en_wflow_G,
      resamples = folds,
      grid = grid_model_G,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with genomic kernel done!\n')
  
  
  
  set.seed(seed)
  xgb_res_GE <-
    tune::tune_grid(
      xgb_wflow_GE,
      resamples = folds,
      grid = xgboost_grid,
      metrics = metric,
      control = tune::control_grid(
        save_pred = TRUE,
        save_workflow = TRUE,
        verbose = FALSE
      )
    )
  cat('Support vector regression with GxE kernel done!\n')
  
  # Initialize a data stack using the stacks() function.
  
  METData_data_st <-
    stacks::stacks() %>%
    stacks::add_candidates(en_res_G) %>%
    stacks::add_candidates(svm_res_E) %>%
    stacks::add_candidates(xgb_res_GE)
  
  # Fit the stack
  
  METData_model_st <-
    METData_data_st %>%
    stacks::blend_predictions()
  
  METData_model_st <-
    METData_model_st %>%
    stacks::fit_members()
  
  # Identify which model configurations were assigned what stacking coefficients
  parameters_collection_G <-
    METData_model_st %>% stacks::collect_parameters('en_res_G')
  parameters_collection_E <-
    METData_model_st %>% stacks::collect_parameters('svm_res_E')
  parameters_collection_GE <-
    METData_model_st %>% stacks::collect_parameters('xgb_res_GE')
  
  
  # Predictions and metrics calculated on a per-environment basis
  
  predictions_test <-
    as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
  
  
  
  # Apply the trained data recipe
  rec_GE <- recipes::prep(rec_GE)
  train_GE = recipes::bake(rec_GE, training)
  test_GE = recipes::bake(rec_GE, test)
  
  
  
  
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'training' = as.data.frame(training),
      'test' = as.data.frame(test)
    ),
    class = c('res_fitted_split','fitted_stacking_reg_3')
  )
  
  
  if (save_model) {
    res_fitted_split[["fitted_model"]] = METData_model_st
  }  else{
    res_fitted_split["fitted_model"] = list(NULL)
  }
  
  
  return(res_fitted_split)
  
  
  
  
}

