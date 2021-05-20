#' @description 
#' @title fit_cv_split
#' bblabbla
#' @param object blabla
#' @export
fit_cv_split <- function(object, ...) {
  UseMethod("fit_cv_split")
}


#' @rdname fit_cv_split
#' @export
fit_cv_split.default <- function(x, ...) {
  stop('not implemented')
  
}

#' @rdname fit_cv_split
#' @export
fit_cv_split.xgb_reg <- function(object,
                                 seed,
                                 inner_cv_reps = 2,
                                 inner_cv_folds = 5,
                                 ...) {
  if (class(object) != "xgb_reg") {
    stop("The object must be an object of the class 'xgb_reg'")
  }
  
  training = object[[1]]
  test = object[[2]]
  rec = object[[3]]
  
  
  
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
  
  # Apply the trained data recipe
  rec <- prep(rec)
  train = bake(rec, training)
  test = bake(rec, test)
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'best_hyperparameters' = best_params,
      'training' = train,
      'test' = test
    ),
    class = 'res_fitted_split'
  )
  
  
  
  
  return(res_fitted_split)
  
  
  
}




#' @rdname fit_cv_split
#' @export
fit_cv_split.xgb_ordinal <- function(object,
                                     seed,
                                     inner_cv_reps = 2,
                                     inner_cv_folds = 5,
                                     ..) {
  if (class(object) != "xgb_ordinal") {
    stop("The object must be an object of the class 'xgb_ordinal'")
  }
  
  training = object[[1]]
  test = object[[2]]
  rec = object[[3]]
  nb_ordinal_classes = object[[4]]
  
  
  # Define the prediction model to use
  
  xgboost_model <-
    parsnip::boost_tree(
      mode = "classification",
      trees = tune(),
      tree_depth = tune(),
      learn_rate = tune()
    ) %>%
    set_engine("xgboost", objective = "binary:logistic") %>%
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
    add_formula(trait_binary ~ .)
  
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
        iter = 8,
        initial = 8,
        metrics = yardstick::metric_set(recall,
                                        precision,
                                        f_meas,
                                        accuracy,
                                        kap,
                                        roc_auc,
                                        sens,
                                        spec),
        control = tune::control_bayes(verbose = FALSE, no_improve = 6)
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
    
  
  
  res_all_binary_classifiers <-
    lapply(seq_len(nb_ordinal_classes - 1), function(x) {
      binary_class_k(
        k = x,
        training = training,
        test = test,
        wf = wf
      )
    })
  
  # Compute the probability of each class based on the results of the serie
  # of binary classifiers
  res_probabilities <-
    as.data.frame(do.call("cbind", res_all_binary_classifiers))
  res_probabilities <-
    res_probabilities %>% select(-contains(c("inf_or_equal")))
  
  list_class_probabilities <- list()
  for (k in seq_len(nb_ordinal_classes)) {
    if (k == 1) {
      list_class_probabilities[[k]] <- 1 - res_probabilities[, 1]
    } else if (k == nb_ordinal_classes) {
      list_class_probabilities[[k]] <- res_probabilities[, k - 1]
    } else {
      list_class_probabilities[[k]] <-
        res_probabilities[, k - 1] - res_probabilities[, k]
    }
  }
  
  
  res_class_probabilities <-
    as.data.frame(do.call("cbind", list_class_probabilities))
  colnames(res_class_probabilities) <-
    paste0('Prob_class=', seq_len(nb_ordinal_classes))
  
  # Obtain the predicted class based on the class probabilities
  
  res_class_probabilities$predicted_class <-
    apply(res_class_probabilities, 1, which.max)
  
  res_class_probabilities <-
    res_class_probabilities %>% bind_cols(test[, trait])
  colnames(res_class_probabilities)[ncol(res_class_probabilities)] <-
    trait
  
  # Compute a serie of performance metrics based on the results from the
  # classification on ordinal data for the test set
  
  res_class_probabilities$predicted_class <-
    as.factor(res_class_probabilities$predicted_class)
  res_class_probabilities[, trait] <-
    as.factor(res_class_probabilities[, trait])
  
  confusion_matrix <-
    yardstick::conf_mat(res_class_probabilities, predicted_class, trait)
  acc <-
    yardstick::accuracy(res_class_probabilities, predicted_class, estimate =
                          get(trait))
  kappa_coefficient <-
    yardstick::kap(res_class_probabilities, predicted_class, estimate = get(trait))
  sens <-
    yardstick::sensitivity(res_class_probabilities, predicted_class, estimate =
                             get(trait))
  spec <-
    yardstick::specificity(res_class_probabilities, predicted_class, estimate =
                             get(trait))
  prec <-
    yardstick::precision(res_class_probabilities, predicted_class, estimate =
                           get(trait))
  recall <-
    yardstick::recall(res_class_probabilities, predicted_class, estimate =
                        get(trait))
  
  # Apply the trained data recipe
  rec <- prep(rec)
  train = bake(rec, training)
  test = bake(rec, test)
  
  # Return final list of class res_fitted_split
  
  res_fitted_split <- structure(
    list(
      'predictions_df' = res_class_probabilities,
      'confusion_matrix' = confusion_matrix,
      'accuracy' = acc,
      'kappa' = kappa_coefficient,
      'sensitivity' = sens,
      'specificity' = spec,
      'precision' = prec,
      'recall' = recall,
      'training' = train,
      'test' = test
    ),
    class = 'res_fitted_split'
  )
  
  
  return(res_fitted_split)
  
  
  
  
  
}



#' @rdname fit_cv_split
#' @export
fit_cv_split.svm_stacking_reg <- function (object,
                                           seed,
                                           inner_cv_reps = 2,
                                           inner_cv_folds = 5,
                                           kernel_G = 'rbf',
                                           kernel_E = 'rbf',
                                           kernel_GE = 'rbf',
                                           ...) {
  # Case if the GE kernel was built (length = 5)
  print(kernel_G)
  print(kernel_E)
  print(kernel_GE)
  
  if (length(object) == 5) {
    training = object[['training']]
    test = object[['test']]
    
    rec_G = object[['rec_G']]
    rec_E = object[['rec_E']]
    rec_GE = object[['rec_GE']]
    
    # Some settings common for all kernels to be trained
    
    metric <- yardstick::metric_set(rmse)
    
    ctrl_res <- stacks::control_stack_resamples()
    
    # Inner CV
    
    set.seed(seed)
    folds <-
      rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    # Define the prediction model to use: different types of kernel can be
    # defined for each of the kernels.
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
    
    #svm_spec_linear <-
    #  parsnip::svm_linear(cost = tune("cost")) %>%
    #  parsnip::set_engine("LiblineaR") %>%
    #  parsnip::set_mode("regression")
    
    
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
    
    if (kernel_GE == 'rbf') {
      svm_spec_GE <- svm_spec_rbf
      grid_model_GE <- 8
    } else if (kernel_GE == 'polynomial') {
      svm_spec_GE <- svm_spec_polynomial
      grid_model_GE <- 14
    } else{
      svm_spec_GE <- svm_spec_linear
      grid_model_GE <- 6
    }
    
    
    # Add recipe and model definition to a workflow, for each of the kernel
    
    svm_wflow_G <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_G) %>%
      workflows::add_recipe(rec_G)
    
    svm_wflow_E <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_E) %>%
      workflows::add_recipe(rec_E)
    
    svm_wflow_GE <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_GE) %>%
      workflows::add_recipe(rec_GE)
    
    # tune cost and sigma with the inner CV for each of the kernels
    
    set.seed(seed)
    svm_res_E <-
      tune::tune_grid(
        svm_wflow_E,
        resamples = folds,
        grid = grid_model_E,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    cat('Support vector regression with env. kernel done!')
    
    set.seed(seed)
    svm_res_GE <-
      tune::tune_grid(
        svm_wflow_GE,
        resamples = folds,
        grid = grid_model_GE,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    cat('Support vector regression with GE kernel done!')
    
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
    cat('Support vector regression with G kernel done!')
    
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
    
    # Identify which model configurations were assigned what stacking coefficients
    parameters_collection_G <-
      METData_model_st %>% stacks::collect_parameters('svm_res_G')
    parameters_collection_E <-
      METData_model_st %>% stacks::collect_parameters('svm_res_E')
    parameters_collection_GE <-
      METData_model_st %>% stacks::collect_parameters('svm_res_GE')
    
    # Predictions
    
    predictions_test <-
      as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      sqrt(mean((
        predictions_test[, trait] - predictions_test[, '.pred']
      ) ^ 2))
    
    # Apply the trained data recipe
    rec_G <- prep(rec_G)
    train_G = bake(rec_G, training)
    test_G = bake(rec_G, test)
    
    rec_E <- prep(rec_E)
    train_E = bake(rec_E, training)
    test_E = bake(rec_E, test)
    
    rec_GE <- prep(rec_GE)
    train_GE = bake(rec_GE, training)
    test_GE = bake(rec_GE, test)
    
    # Return final list of class res_fitted_split
    
    res_fitted_split <- structure(
      list(
        'parameters_collection_G' = parameters_collection_G,
        'parameters_collection_E' = parameters_collection_E,
        'parameters_collection_GE' = parameters_collection_GE,
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'training_G' = train_G,
        'training_E' = train_E,
        'training_GE' = train_GE,
        'test_G' = test_G,
        'test_E' = test_E,
        'test_GE' = test_GE
      ),
      class = 'res_fitted_split'
    )
  }
  
  # Case if the GE kernel was not built (length = 4)
  
  if (length(object) == 4) {
    training = object[['training']]
    test = object[['test']]
    
    rec_G = object[['rec_G']]
    rec_E = object[['rec_E']]
    
    
    # Some settings common for all kernels to be trained
    
    metric <- yardstick::metric_set(rmse)
    
    ctrl_res <- stacks::control_stack_resamples()
    
    # Inner CV
    
    set.seed(seed)
    folds <-
      rsample::vfold_cv(training, repeats = inner_cv_reps, v = inner_cv_folds)
    
    # Define the prediction model to use: different types of kernel can be
    # defined for each of the kernels.
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
    
    #svm_spec_linear <-
    #  parsnip::svm_linear(cost = tune("cost")) %>%
    #  parsnip::set_engine("LiblineaR") %>%
    #  parsnip::set_mode("regression")
    
    
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
    
    if (kernel_GE == 'rbf') {
      svm_spec_GE <- svm_spec_rbf
      grid_model_GE <- 8
    } else if (kernel_GE == 'polynomial') {
      svm_spec_GE <- svm_spec_polynomial
      grid_model_GE <- 14
    } else{
      svm_spec_GE <- svm_spec_linear
      grid_model_GE <- 6
    }
    
    
    # Add recipe and model definition to a workflow, for each of the kernel
    
    svm_wflow_G <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_G) %>%
      workflows::add_recipe(rec_G)
    
    svm_wflow_E <-
      workflows::workflow() %>%
      workflows::add_model(svm_spec_E) %>%
      workflows::add_recipe(rec_E)
    
    
    # tune cost and sigma with the inner CV for each of the kernels
    
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
    
    # Initialize a data stack using the stacks() function.
    
    METData_data_st <-
      stacks::stacks() %>%
      stacks::add_candidates(svm_res_G) %>%
      stacks::add_candidates(svm_res_E)
    
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
    
    # Predictions
    
    predictions_test <-
      as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
    rmse_pred_obs <-
      sqrt(mean((
        predictions_test[, trait] - predictions_test[, '.pred']
      ) ^ 2))
    
    # Apply the trained data recipe
    rec_G <- prep(rec_G)
    train_G = bake(rec_G, training)
    test_G = bake(rec_G, test)
    
    rec_E <- prep(rec_E)
    train_E = bake(rec_E, training)
    test_E = bake(rec_E, test)
    
    
    # Return final list of class res_fitted_split
    
    res_fitted_split <- structure(
      list(
        'parameters_collection_G' = parameters_collection_G,
        'parameters_collection_E' = parameters_collection_E,
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'training_G' = train_G,
        'training_E' = train_E,
        'test_G' = test_G,
        'test_E' = test_E
      ),
      class = 'res_fitted_split'
    )
  }
  
  return(res_fitted_split)
  
  
}




#' @rdname fit_cv_split
#' @export
fit_cv_split.xgb_multiclass_factor <- function(object,
                                               seed,
                                               inner_cv_reps = 2,
                                               inner_cv_folds = 5,
                                               ..) {
  if (class(object) != "xgb_multiclass_factor") {
    stop("The object must be an object of the class 'xgb_multiclass_factor'")
  }
  
  training = object[[1]]
  test = object[[2]]
  rec = object[[3]]
  nb_classes = object[[4]]
  
  
  # Define the prediction model to use
  
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
    add_formula(trait ~ .)
  
  
    
    
    
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
        iter = 8,
        initial = 8,
        metrics = yardstick::metric_set(recall,
                                        precision,
                                        f_meas,
                                        accuracy,
                                        kap,
                                        roc_auc,
                                        sens,
                                        spec),
        control = tune::control_bayes(verbose = FALSE, no_improve = 6)
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
    lapply(seq_len(nb_ordinal_classes - 1), function(x) {
      binary_class_k(
        k = x,
        training = training,
        test = test,
        wf = wf
      )
    })
  
  # Compute the probability of each class based on the results of the serie
  # of binary classifiers
  res_probabilities <-
    as.data.frame(do.call("cbind", res_all_binary_classifiers))
  res_probabilities <-
    res_probabilities %>% select(-contains(c("inf_or_equal")))
  
  list_class_probabilities <- list()
  for (k in seq_len(nb_ordinal_classes)) {
    if (k == 1) {
      list_class_probabilities[[k]] <- 1 - res_probabilities[, 1]
    } else if (k == nb_ordinal_classes) {
      list_class_probabilities[[k]] <- res_probabilities[, k - 1]
    } else {
      list_class_probabilities[[k]] <-
        res_probabilities[, k - 1] - res_probabilities[, k]
    }
  }
  
  
  res_class_probabilities <-
    as.data.frame(do.call("cbind", list_class_probabilities))
  colnames(res_class_probabilities) <-
    paste0('Prob_class=', seq_len(nb_ordinal_classes))
  
  # Obtain the predicted class based on the class probabilities
  
  res_class_probabilities$predicted_class <-
    apply(res_class_probabilities, 1, which.max)
  
  res_class_probabilities <-
    res_class_probabilities %>% bind_cols(test[, trait])
  colnames(res_class_probabilities)[ncol(res_class_probabilities)] <-
    trait
  
  # Compute a serie of performance metrics based on the results from the
  # classification on ordinal data for the test set
  
  res_class_probabilities$predicted_class <-
    as.factor(res_class_probabilities$predicted_class)
  res_class_probabilities[, trait] <-
    as.factor(res_class_probabilities[, trait])
  
  confusion_matrix <-
    yardstick::conf_mat(res_class_probabilities, predicted_class, trait)
  acc <-
    yardstick::accuracy(res_class_probabilities, predicted_class, estimate =
                          get(trait))
  kappa_coefficient <-
    yardstick::kap(res_class_probabilities, predicted_class, estimate = get(trait))
  sens <-
    yardstick::sensitivity(res_class_probabilities, predicted_class, estimate =
                             get(trait))
  spec <-
    yardstick::specificity(res_class_probabilities, predicted_class, estimate =
                             get(trait))
  prec <-
    yardstick::precision(res_class_probabilities, predicted_class, estimate =
                           get(trait))
  recall <-
    yardstick::recall(res_class_probabilities, predicted_class, estimate =
                        get(trait))
  
  # Apply the trained data recipe
  rec <- prep(rec)
  train = bake(rec, training)
  test = bake(rec, test)
  
  # Return final list of class res_fitted_split
  
  res_fitted_split <- structure(
    list(
      'predictions_df' = res_class_probabilities,
      'confusion_matrix' = confusion_matrix,
      'accuracy' = acc,
      'kappa' = kappa_coefficient,
      'sensitivity' = sens,
      'specificity' = spec,
      'precision' = prec,
      'recall' = recall,
      'training' = train,
      'test' = test
    ),
    class = 'res_fitted_split'
  )
  
  
  return(res_fitted_split)
  
  
  
  
  
}

# == 'naive_multiclass') {
#  xgboost_model <-
#    parsnip::boost_tree(
#      mode = "classification",
#      trees = tune(),
#      tree_depth = tune(),
#      learn_rate = tune()
#    ) %>%
#    set_engine(
#      "xgboost",
#      objective = "multi:softprob",
#      lambda = 0,
#      alpha = 1,
#      num_class = nb_classes
#    ) %>%
#    translate()
#} 