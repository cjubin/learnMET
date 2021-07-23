

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

