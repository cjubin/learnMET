#' @rdname fit_cv_split
#' @export
fit_cv_split.DL_reg <- function(object,
                                seed,
                                inner_cv_reps = 1,
                                inner_cv_folds = 3,
                                vip = F,
                                ...) {
  if (class(object) != "DL_reg") {
    stop("The object must be an object of the class 'DL_reg'")
  }

  training = object[['training']]
  test = object[['test']]
  test_before_recipe = test
  rec = object[['rec']]
  trait = as.character(rec$var_info[rec$var_info$role == 'outcome', 'variable'])
  
  prepped_recipe <- prep(rec, strings_as_factors = FALSE)
  training <- bake(prepped_recipe, training)
  test <- bake(prepped_recipe, test)
  
  all_predictors <-
    as.character(prepped_recipe$term_info[prepped_recipe$term_info$role == 'predictor', 'variable']$variable)
  
  ## Split of the training set in a training and validation set to optimize
  ## hyperparameters
  
  split_tr <- initial_split(training, prop = 0.6)
  split_tr$out_id <-
    (1:nrow(training))[(1:nrow(training)) %notin% split_tr$in_id]
  # training set
  tr_data_x <-
    as.matrix((training %>%
                 dplyr::select(-IDenv,-all_of(trait)))[split_tr$in_id,])
  tr_data_y <-
    as.matrix((training %>%
                 dplyr::select(all_of(trait)))[split_tr$in_id,])
  # validation set
  val_data_x <-
    as.matrix((training %>%
                 dplyr::select(-IDenv,-all_of(trait)))[split_tr$out_id,])
  val_data_y <-
    as.matrix((training %>%
                 dplyr::select(all_of(trait)))[split_tr$out_id,])
  
  # Define the prediction model to use
  
  keras_fit <-
    function(units_1,
             units_2,
             dropout1,
             dropout2,
             learning_rate) {
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
          optimizer = optimizer_adam(lr = learning_rate),
          metrics = list("mean_absolute_error")
        )
      
      history <- DL_model %>% fit(
        tr_data_x,
        tr_data_y,
        batch_size = 128,
        epochs = 250,
        verbose = 0,
        validation_data = list(val_data_x, val_data_y)
      )
      
      print(names(history$metrics))
      
      result <-
        list(Score = -history$metrics$val_mean_absolute_error[250],
             Pred = 0)
      
      return(result)
      
    }
  search_bound_keras <- list(
    units_1 = c(40, 80),
    units_2 = c(20, 60),
    dropout1 = c(0.1, 0.5),
    dropout2 = c(0.1, 0.5),
    learning_rate = c(0.001, 0.01)
  )
  
  
  set.seed(seed)
  bayes_keras <-
    rBayesianOptimization::BayesianOptimization(
      FUN = keras_fit,
      bounds = search_bound_keras,
      init_points = 15,
      init_grid_dt = NULL,
      n_iter = 20,
      acq = "ucb"
    )
  
  
  
  
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
    layer_dense(units = ceiling(bayes_keras$Best_Par['units_2']),
                activation = 'relu') %>%
    layer_dropout(rate = bayes_keras$Best_Par['dropout2']) %>%
    layer_dense(units = 1, activation = 'linear') %>%
    compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(lr = bayes_keras$Best_Par['learning_rate']),
      metrics = list("mean_absolute_error")
    )
  
  best_params = bayes_keras$Best_Par
  
  
  DL_model %>% fit(
    x = as.matrix(training %>%
                    dplyr::select(-IDenv,-all_of(trait))),
    y = as.matrix(training %>%
                    dplyr::select(all_of(trait))),
    batch_size = 128,
    epochs = 250,
    verbose = 0,
    validation_split = 0.3
  )
  
  predictions_test <-
    as.data.frame(DL_model %>%
                    predict(x = as.matrix(
                      test %>%
                        dplyr::select(-IDenv,-all_of(trait))
                    )))
  colnames(predictions_test) <- '.pred'
  predictions_test <-
    predictions_test %>% bind_cols(test_before_recipe)
  
  cor_pred_obs <-
    predictions_test %>%
    group_by(IDenv) %>% summarize(COR = cor(.pred, get(trait), method = 'pearson'))
  
  rmse_pred_obs <-
    predictions_test %>%
    group_by(IDenv) %>% summarize(RMSE = sqrt(mean((get(
      trait
    ) - .pred) ^ 2)))
  
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'best_hyperparameters' = as.data.frame(best_params),
      'training' = as.data.frame(training),
      'test' = as.data.frame(test),
      'vip' = data.frame()
    ),
    class = 'res_fitted_split'
  )
  
   
  if (vip){
  fitted_obj_for_vip <- structure(
    list(
      model = DL_model,
      x_train = as.matrix(training %>%
                            dplyr::select(-IDenv,-all_of(trait))),
      y_train = as.matrix(training %>%
                            dplyr::select(all_of(trait))),
      trait = trait
    ),
    class = c('fitted_DL_reg', 'list')
  )
  
  # Obtain the variable importance
  
  variable_importance_vip <-
    variable_importance_split(fitted_obj_for_vip)
  
  
  
  # Return final list of class res_fitted_split
  res_fitted_split <- structure(
    list(
      'prediction_method' = class(object),
      'predictions_df' = predictions_test,
      'cor_pred_obs' = cor_pred_obs,
      'rmse_pred_obs' = rmse_pred_obs,
      'best_hyperparameters' = as.data.frame(best_params),
      'training' = as.data.frame(training),
      'test' = as.data.frame(test),
      'vip' = variable_importance_vip
    ),
    class = 'res_fitted_split'
  )
  
  
}
  
  
  
  
  
  return(res_fitted_split)
  
  
  
}
