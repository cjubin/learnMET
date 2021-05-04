#' Compute marker effects per environment with Elastic Net
#'
#' @param environment \code{character} indicating the name of the environment
#'   for which marker effects should be computed
#'
#' @param geno \code{data.frame} with markers in columns and inviduals in rows.
#'   Typical input is METData$geno.
#'
#' @param pheno \code{data.frame} with:
#'   First column: ID genotypes
#'   Second column: year
#'   Third column: location
#'   Subsequent columns: phenotypic traits with names indicated in colnames()
#'   Last column: IDenv (combination LocationxYear)
#'   Typical input is METData$pheno
#'
#' @param pheno_trait \code{character} Name of the trait under study for which 
#' marker effects should be estimated
#'
#' @param nb_folds_cv \code{numeric} Number of folds in the CV process
#'
#' @param reps \code{numeric} Number of repeats of the k-folds CV
#'
#' @return all_coef_avg \code{data.frame}
#'   First column \code{character} contains the marker names.
#'   Second column \code{numeric} the marker effects in this environment
#'   calculated by cross-validation.
#'   Third column \code{character} contains the environment name (combination
#'   LocationxYear).
#'


marker_effect_per_env_EN <-
  
  function(geno,
           pheno,
           environment,
           pheno_trait,
           nb_folds_cv = 5,
           reps = 2,
           ...) {
    # Select the phenotype data corresponding to the selected environment
    
    pheno <- pheno[pheno$IDenv == environment, ]
    list_predictors <- colnames(geno)
    geno$geno_ID = row.names(geno)
    pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)
    
    # Select trait and marker columns from the phenotypic file merged with geno
    # data
    
    pheno <- pheno[, c(pheno_trait, list_predictors)]
    
    #Create the cross-validation random splits
    
    cv_splits <-
      rsample::vfold_cv(pheno, folds = nb_folds_cv, repeats = reps)
    
    # Define the type of model to tune
    
    mod <- parsnip::linear_reg(penalty = tune(),
                               mixture = tune()) %>%
      parsnip::set_engine("glmnet")
    
    
    # Define the predictors and outcome variables.
    # Remove predictors with null varaiance.
    # Center and scale the variables (standardization)
    
    rec <- recipes::recipe(~ ., data = pheno) %>%
      recipes::update_role(all_of(pheno_trait), new_role = 'outcome') %>%
      recipes::update_role(all_of(list_predictors), new_role = 'predictor') %>%
      recipes::step_nzv(all_numeric()) %>%
      recipes::step_normalize(all_numeric())
    
    # Define a workflow based on an Elastic net model
    
    wfl <- workflows::workflow() %>%
      workflows::add_recipe(rec) %>%
      workflows::add_model(mod)
    
    mixture_param <-
      tune::parameters(dials::penalty(), dials::mixture())
    
    # Define the grid search space for tuning parameters with the range of values
    # for penalty and mixture coefficient
    
    mixture_param  <-
      mixture_param  %>% update(penalty = dials::penalty(range = c(-5, 1), trans = scales::log10_trans()))
    
    glmn_grid <-
      dials::grid_regular(mixture_param, levels = c(10, 10))
    
    
    ctrl <- tune::control_grid(save_pred = TRUE, verbose = F)
    
    # Tuning hyperparameters with the defined resampling strategy and grid of
    # hyperparameters
    
    glmn_tune <-
      tune::tune_grid(
        wfl,
        resamples = cv_splits,
        grid = glmn_grid,
        metrics = yardstick::metric_set(rmse),
        control = ctrl
      )
    
    # Select the best hyperparameters for estimating marker effects
    # by cross-validation
    
    best_parameters <-
      tune::select_best(glmn_tune, metric = "rmse")
    
    wfl_final <-
      wfl %>%
      tune::finalize_workflow(best_parameters)
    
    
    get_model <- function(x) {
      pull_workflow_fit(x) %>% tidy()
    }
    ctrl <- control_resamples(extract = get_model)
    
    glmn_cv_final <-
      tune::fit_resamples(wfl_final, cv_splits, control = ctrl)
    
    all_coef <- map_dfr(glmn_cv_final$.extracts, ~ .x[[1]][[1]])
    
    # Note: markers for which the estimated effect is 0 are not included in this
    # table. To calculate the average effect, we need to divide the sum of marker
    # estimates per marker by the number of reps*nb_folds_cv.
    
    all_coef_avg <-
      all_coef %>% group_by(term) %>% summarise(cv_mean = sum(estimate) / (reps *
                                                                             nb_folds_cv))
    
    #all_coef_avg$abs_cv_mean <- abs(all_coef_avg$cv_mean)
    all_coef_avg$environment <- environment
    
    
    return(all_coef_avg)
    
  }
