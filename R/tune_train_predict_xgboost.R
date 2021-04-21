# Tuning hyperparameters for each train/test split with optimization on inner CV

tune_train_predict_xgboost <-
  function(training_set,
           test_set,
           inner_nb_folds = 2,
           inner_reps = 1,
           bayesian_optimization = T,
           grid_search_optimization = F,
           pheno_trait = trait,
           ...) {
    # Create the cross-validation random splits for hyperparameter optimization
    
    cv_splits <-
      rsample::vfold_cv(training_set, folds = inner_nb_folds, repeats = inner_reps)
    
    # Define the type of model to tune
    
    xgboost_model <- parsnip::boost_tree(
      mode = "regression",
      trees = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      sample_size = 1
    ) %>%
      set_engine("xgboost", objective = "reg:linear") %>%
      translate()
    
    
    xgb_grid <- parameters(trees(),
                           learn_rate(),
                           tree_depth())
    
    # Define the formula based on the trait to study
    # All predictors previously defined to be used are included.
    
    form <- as.formula(paste(trait, ' ~ .'))
    
    wf <- workflow() %>%
      add_model(xgboost_model) %>%
      add_formula(formula = form)
    
    if (bayesian_optimization) {
      opt_res <- wf %>%
        tune_bayes(
          resamples = cv_splits,
          param_info = xgb_grid,
          iter = 25,
          initial = 10,
          metrics = yardstick::metric_set(rmse),
          control = tune::control_bayes(verbose = TRUE, no_improve = 15)
        )
      
    }
    
    #The training-inner data set was used to train the DL model using the grid of hyperparameters values. This inner CV strategy was facilitated by using the internal capabilities of Keras by means of the validation_split argument on the fit() function. The predictve power is assessed in the second part of the data set (testing-inner). With this, a set of best-fitting hyperparameters (best combination of units, epoch and layers) from the inner CV loop is obtained. Finally, these set of hyperparameters were used to predict the performance in the independent testing data set (testing-outer).
    
    # Keep the best hyperparameters for model training
    
    xgboost_best_params <- opt_res %>%
      tune::select_best("rmse", maximize = FALSE)
    
    
    xgboost_model_final <- finalize_workflow(wf,
                                             xgboost_best_params)
    
    final_xgboost <-
      xgboost_model_final %>%
      fit(data = training_set)
    
    
  }
