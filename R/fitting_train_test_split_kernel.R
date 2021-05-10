#' Multiple kernel learning method based on G, E and GxE datasets.
#'
#' @description
#' 
#' A multiple kernel support vector machine framework using model stacking.
#' Model stacking is an ensemble method that takes the outputs of different 
#' support vector machine models constructed based on different data types 
#' (genomic, environmental and GxE interactions). For each type of kernel, 
#' multiple model configurations are defined based on a grid of hyperparameters.
#' Then, an ensemble is built with stacks to create an object that contain the 
#' assessment set predictions for each candidate ensemble member. A LASSO model 
#' is used to figure out how the respective output of each model from the stack 
#' members should be combined to obtain a final prediction, and to estimate the 
#' "stacking coefficients" of the model stack. Candidate members with a stacking 
#' coefficient different from 0 are trained on the full training set, and the 
#' test set can be predicted using the "instructions" on how to combine the
#' respective predictions.
#' 
#' @details For more information, consult: 
#'  \url{https://stacks.tidymodels.org/index.html}
#' 
#'
#' @param split. A \code{split_processed} object containing:
#' \describe{
#'   \item{training}{\code{data.frame} Training set.}
#'   \item{test}{\code{data.frame} Test set.}
#'   \item{rec_G}{\code{recipe} object with the different steps to process 
#'    molecular marker dataset from the training set. Same transformations 
#'    applied on the test set.}
#'   \item{rec_E}{\code{recipe} object with the different steps to process 
#'    environmental data from the training set. Same transformations applied on
#'    the test set.}
#'   \item{rec_GE}{\code{recipe} object with the different steps to process 
#'    the GxE dataset from the training set. Same transformations applied on
#'    the test set.}
#'    }
#' @param seed \code{integer} Seed value.
#' 
#' @param inner_cv_reps \code{integer} Number of times to repeat the k-fold 
#'   partitioning used for the inner cross-validation for estimation of the best 
#'   hyperparameters. The same resampling object is used for all kernel 
#'   configurations which are evaluated. Default is 2.
#' 
#' @param inner_cv_folds \code{integer} Number of partitions of the training set
#'   for the inner cross-validation for estimation of the best hyperparameters.
#'   The same resampling object is used for all kernel 
#'   configurations which are evaluated. Default is 5.
#'   
#' @param kernel_G \code{character} Type of kernel function to use for the 
#' molecular marker dataset. Options are `rbf` (default), `polynomial` or 
#' `linear`.
#'  
#' @param kernel_E \code{character} Type of kernel function to use for the
#' environmental dataset. Options are `rbf` (default), `polynomial` or `linear`.
#'   
#' @param kernel_GE \code{character} Type of kernel function to use for the GxE
#' dataset. Options are `rbf` (default), `polynomial` or `linear`.
#' 
#' @return a \code{list} object containing:
#'   \describe{
#'     \item{training}{\code{data.frame} Training set.}
#'     \item{test}{\code{data.frame} Test set.}
#'     \item{parameters_collection_G}{a \code{tbl_df} mapping all of the 
#'      candidate models based on genomic data to their hyperparameters with 
#'      their stacking coefficient obtained after evaluation of the data stack 
#'      on the training set.}
#'     \item{parameters_collection_E}{a \code{tbl_df} mapping all of the 
#'      candidate models based on environmental data to their hyperparameters 
#'      with their stacking coefficient.}
#'     \item{parameters_collection_GE}{a \code{tbl_df} mapping all of the 
#'      candidate models based on GxE data to their hyperparameters with their 
#'      stacking coefficient.}
#'     \item{predictions_df}{\code{data.frame} with original test dataset with 
#'      extra column containing predicted values.}
#'     \item{cor_pred_obs}{\code{numeric} Pearson's correlation between 
#'      predicted and observed values of the test set.}
#'     \item{rmse_pred_obs}{\code{numeric} root mean square error between 
#'      predicted and observed values of the test set.}
#'  }
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
#'



fitting_train_test_split_kernel <-
  function(split,
           seed,
           inner_cv_reps = 2,
           inner_cv_folds = 5,
           kernel_G = 'rbf',
           kernel_E = 'rbf',
           kernel_GE = 'rbf',
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
        grid = grid_model_G,
        metrics = metric,
        control = tune::control_grid(save_pred = TRUE,
                                     save_workflow = TRUE)
      )
    
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
      sqrt(mean((predictions_test[, trait] - predictions_test[, '.pred']) ^ 2))
    
    
    return
    
    
    return(
      list(
        'training' = training,
        'test' = test,
        'parameters_collection_G' = parameters_collection_G,
        'parameters_collection_E' = parameters_collection_E,
        'parameters_collection_GE' = parameters_collection_GE,
        'predictions_df' = predictions_test,
        'cor_pred_obs' = cor_pred_obs,
        'rmse_pred_obs' = rmse_pred_obs
        
      )
    )
    
    
    
  }