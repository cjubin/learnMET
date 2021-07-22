#' @description
#' for all CV schemes, the prediction accuracy is computed as the correlations between the observed
#' and predicted values within same environments.
#' @title fit_cv_split
#'
#' @param fit CV object
#' @param object It must be an object of class
#'
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
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
                                 inner_cv_reps = 1,
                                 inner_cv_folds = 3,
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
      iter = 5,
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
  
  # Obtain the variable importance with Shapley importance
  
  xg_mod <- model_final %>%
    fit(data = training) %>%
    pull_workflow_fit()
  
  X <- prep(rec, training) %>%
    juice() %>%
    dplyr::select(-IDenv) %>%
    dplyr::select(-all_of(trait)) %>%
    as.data.frame() %>%
    as.matrix()
  
  shap <- fastshap::explain(xg_mod$fit, X = X, exact = TRUE)
  shap_data <- autoplot(shap)$data
  
  
  
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
      'ranking_vip' = ranking_vip,
      'shapley_importance' = shap_data
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
                                     ...) {
  if (class(object) != "xgb_ordinal") {
    stop("The object must be an object of the class 'xgb_ordinal'")
  }
  
  training = object[['training']]
  test = object[['test']]
  rec = object[['rec']]
  trait = as.character(rec$var_info[rec$var_info$role == 'outcome', 'variable'])
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
    
    prepped <- prep(rec,strings_as_factors = FALSE)
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
    res_class_probabilities %>% bind_cols(test[, trait]) %>% bind_cols(test[, 'IDenv'])
  
  colnames(res_class_probabilities)[ncol(res_class_probabilities)] <-
    trait
  
  
  
  # Compute a serie of performance metrics based on the results from the
  # classification on ordinal data for the test set
  
  res_class_probabilities$predicted_class <-
    as.factor(res_class_probabilities$predicted_class)
  res_class_probabilities[, trait] <-
    as.factor(res_class_probabilities[, trait])
  
  confusion_matrix <-
    res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::conf_mat(estimate = predicted_class, truth = trait)
  
  acc <- res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::accuracy(estimate = predicted_class,
                        truth =
                          get(trait))
  
  kappa_coefficient <-
    res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::kap(estimate = predicted_class,
                   truth = get(trait))
  
  sens <- res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::sensitivity(estimate = predicted_class,
                           truth = get(trait))
  
  spec <- res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::specificity(estimate = predicted_class,
                           truth = get(trait))
  
  prec <- res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::precision(estimate = predicted_class,
                         truth = get(trait))
  recall <- res_class_probabilities %>% group_by(IDenv) %>%
    yardstick::recall(estimate = predicted_class,
                      truth = get(trait))
  
  # Apply the trained data recipe
  rec <- prep(rec,strings_as_factors = FALSE)
  train = bake(rec, training)
  test = bake(rec, test)
  
  # Return final list of class res_fitted_split
  
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = res_class_probabilities,
      'confusion_matrix' = confusion_matrix,
      'accuracy' = acc,
      'kappa' = kappa_coefficient,
      'sensitivity' = sens,
      'specificity' = spec,
      'precision' = prec,
      'recall' = recall,
      'training_transformed' = as.data.frame(train),
      'test_transformed' = as.data.frame(test)
    ),
    class = 'res_fitted_split'
  )
  
  
  return(res_fitted_split)
  
  
}


#' @rdname fit_cv_split
#' @export
fit_cv_split.DL_reg <- function(object,
                                seed,
                                inner_cv_reps = 1,
                                inner_cv_folds = 3,
                                ...) {
  if (class(object) != "DL_reg") {
    stop("The object must be an object of the class 'DL_reg'")
  }
  
  training = object[['training']]
  test = object[['test']]
  test_before_recipe = test
  rec = object[['rec']]
  trait = as.character(rec$var_info[rec$var_info$role == 'outcome', 'variable'])
  
  prepped_recipe <- prep(rec,strings_as_factors = FALSE)
  training <- bake(prepped_recipe,training)
  test <- bake(prepped_recipe,test)
  
  all_predictors <-
    as.character(prepped_recipe$term_info[prepped_recipe$term_info$role == 'predictor', 'variable']$variable)
  
  ## Split of the training set in a training and validation set to optimize
  ## hyperparameters

  split_tr <- initial_split(training,prop = 0.6)
  split_tr$out_id <- (1:nrow(training))[(1:nrow(training))%notin%split_tr$in_id]
  # training set
  tr_data_x <- as.matrix((training %>% select(-IDenv,-all_of(trait)))[split_tr$in_id,])
  tr_data_y <- as.matrix((training %>% select(all_of(trait)))[split_tr$in_id,])
  # validation set
  val_data_x <- as.matrix((training %>% select(-IDenv,-all_of(trait)))[split_tr$out_id,])
  val_data_y <- as.matrix((training %>% select(all_of(trait)))[split_tr$out_id,])
  
  # Define the prediction model to use
  
  keras_fit <- function(units_1,units_2, dropout1,dropout2, learning_rate){
    
  DL_model <- keras_model_sequential() %>%
    layer_dense(
      units = ceiling(units_1),
      activation = 'relu',
      input_shape = c(ncol(tr_data_x))
    ) %>%
    layer_dropout(rate = dropout1) %>%
    layer_dense(units = ceiling(units_2), activation = 'relu') %>%
    layer_dropout(rate = dropout2) %>%
    layer_dense(units = 1, activation = 'linear') %>%
    compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr=learning_rate),
      metrics = list("mean_absolute_error")
    )
  
  history <- DL_model %>% fit(
    tr_data_x, tr_data_y,
    batch_size = 128, 
    epochs = 150,
    verbose = 0,
    validation_data = list(
      val_data_x,val_data_y
    ))
    
    print(names(history$metrics))
    
    result <- list(Score = -history$metrics$val_mean_absolute_error[150], Pred = 0)
                   
    return(result)
  
  }
  search_bound_keras <- list(units_1 = c(40,100),
                             units_2 = c(20,50),
                             dropout1 = c(0,0.5),
                             dropout2 = c(0,0.5),
                             learning_rate = c(0.001, 0.01))
  
 
  set.seed(seed)
  bayes_keras <- rBayesianOptimization::BayesianOptimization(FUN = keras_fit, bounds = search_bound_keras, 
                                      init_points = 15, init_grid_dt = NULL, 
                                      n_iter = 20, acq = "ucb")
  
  
 
  
  cat('Optimizing hyperparameters for this training set: done!\n')
  
  # Retain the best hyperparameters and update the workflow with these
  # hyperparameters
  
  DL_model <- keras_model_sequential() %>%
    layer_dense(
      units = ceiling(bayes_keras$Best_Par['units_1']),
      activation = 'relu',
      input_shape = c(ncol(tr_data_x))
    ) %>%
    layer_dropout(rate = bayes_keras$Best_Par['dropout1']) %>%
    layer_dense(units = ceiling(bayes_keras$Best_Par['units_2']), activation = 'relu') %>%
    layer_dropout(rate = bayes_keras$Best_Par['dropout2']) %>%
    layer_dense(units = 1, activation = 'linear') %>%
    compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr=bayes_keras$Best_Par['learning_rate']),
      metrics = list("mean_absolute_error")
    )
  
   best_params = bayes_keras$Best_Par
  
  
  DL_model %>% fit(
    x = as.matrix(training %>% select(-IDenv,-all_of(trait))),y=as.matrix(training %>% select(all_of(trait))),
    batch_size = 128, 
    epochs = 70,
    verbose = 0,
    validation_split = 0.2)
  
  predictions_test <- as.data.frame(DL_model %>% predict(x = as.matrix(test %>% select(-IDenv,-all_of(trait)))))
  colnames(predictions_test)<- '.pred'
  predictions_test <- predictions_test %>% bind_cols(test_before_recipe)
  
  cor_pred_obs <-
    predictions_test %>%
    group_by(IDenv) %>% summarize(COR = cor(.pred, get(trait), method = 'pearson'))
  
  rmse_pred_obs <-
    predictions_test %>%
    group_by(IDenv) %>% summarize(RMSE = sqrt(mean((get(
      trait
    ) - .pred) ^ 2)))
  
  
  # Obtain the variable importance
  
  pred_wrapper <- function(object, newdata) {
    predict(object, x = as.matrix(newdata)) %>%
      as.vector()
  }
  
  variable_importance_vip <- vip(
    object = DL_model,                     
    method = "permute",                 
    num_features = ncol(tr_data_x),      
    pred_wrapper = pred_wrapper,           
    target = tr_data_y,            
    metric = "rsquared",               
    train = as.data.frame(tr_data_x)
  )
  
   
  ranking_vip <- as.data.frame(variable_importance_vip$data)
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'best_hyperparameters' = as.data.frame(best_params),
      'training_transformed' = as.data.frame(training),
      'test_transformed' = as.data.frame(test),
      'ranking_vip' = ranking_vip
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
    trait = as.character(rec_G$var_info[rec_G$var_info$role == 'outcome', 'variable'])
    
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
      grid_model_G <- 6
    } else if (kernel_G == 'polynomial') {
      svm_spec_G <- svm_spec_polynomial
      grid_model_G <- 14
    } else{
      svm_spec_G <- svm_spec_linear
      grid_model_G <- 6
    }
    
    if (kernel_E == 'rbf') {
      svm_spec_E <- svm_spec_rbf
      grid_model_E <- 6
    } else if (kernel_E == 'polynomial') {
      svm_spec_E <- svm_spec_polynomial
      grid_model_E <- 14
    } else{
      svm_spec_E <- svm_spec_linear
      grid_model_E <- 6
    }
    
    if (kernel_GE == 'rbf') {
      svm_spec_GE <- svm_spec_rbf
      grid_model_GE <- 6
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
    
    # Predictions and metrics calculated on a per-environment basis
    
    predictions_test <-
      as.data.frame(METData_model_st %>% predict(new_data = test) %>% bind_cols(test))
    
    cor_pred_obs <-
      METData_model_st %>% predict(new_data = test) %>% bind_cols(test) %>%
      group_by(IDenv) %>% summarize(COR = cor(.pred, get(trait), method = 'pearson'))
    print(cor_pred_obs)
    
    #cor(predictions_test[, '.pred'], predictions_test[, trait], method = 'pearson')
    
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
        'prediction_method' = class(object),
        'parameters_collection_G' = as.data.frame(parameters_collection_G),
        'parameters_collection_E' = as.data.frame(parameters_collection_E),
        'parameters_collection_GE' = as.data.frame(parameters_collection_GE),
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'training_G_transformed' = as.data.frame(train_G),
        'training_E_transformed' = as.data.frame(train_E),
        'training_GE_transformed' = as.data.frame(train_GE),
        'test_G' = as.data.frame(test_G),
        'test_E' = as.data.frame(test_E),
        'test_GE' = as.data.frame(test_GE)
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
    trait = as.character(rec_G$var_info[rec_G$var_info$role == 'outcome', 'variable'])
    
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
    rec_G <- prep(rec_G)
    train_G = bake(rec_G, training)
    test_G = bake(rec_G, test)
    
    rec_E <- prep(rec_E)
    train_E = bake(rec_E, training)
    test_E = bake(rec_E, test)
    
    
    # Return final list of class res_fitted_split
    
    res_fitted_split <- structure(
      list(
        'prediction_method' = class(object),
        'parameters_collection_G' = as.data.frame(parameters_collection_G),
        'parameters_collection_E' = as.data.frame(parameters_collection_E),
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs,
        'training_G_transformed' = as.data.frame(train_G),
        'training_E_transformed' = as.data.frame(train_E),
        'test_G' = as.data.frame(test_G),
        'test_E' = as.data.frame(test_E)
      ),
      class = 'res_fitted_split'
    )
  }
  
  return(res_fitted_split)
  
  
}
