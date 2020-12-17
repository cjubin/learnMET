#default elasticnet
#but add rqtl + ridge reg ,and RF/xgboost later


select_markers <- function(METData,
                           trait,
                           method = 'elasticnet',
                           nb_folds_cv = 10,
                           reps = 1,
                           nb_selected_markers = 1000) {
  ##

  if (METData$filtering_markers == TRUE &&
      dim(METData$geno[2] < 1000) {
        stop('The number of markers is low and does not need to be further reduced.')
      }



      ## Separate analysis per environment

      marker_effect_per_env <- function(environment, geno, pheno,trait) {

        pheno <- pheno[pheno$IDenv == environment,]
        geno <- as.data.frame(geno)
        list_predictors <- colnames(geno)
        geno$geno_ID = row.names(geno)
        pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)

        # Reduce pheno data.frame to trait +  markers as columns

        pheno <- pheno[, c(trait, list_predictors)]

        ## Create the cross-validation randomly splits

        cv_splits <-
          rsample::vfold_cv(pheno, folds = nb_folds_cv, repeats = reps)

        mod <- parsnip::linear_reg(penalty = tune(),
                                   mixture = tune()) %>%
          set_engine("glmnet")



        rec <- recipe( ~ ., data = pheno) %>%
          update_role(all_of(trait), new_role = 'outcome') %>%
          update_role(all_of(list_predictors), new_role = 'predictor') %>%
          step_nzv(all_numeric()) %>%
          step_normalize(all_numeric())

        wfl <- workflow() %>%
          add_recipe(rec) %>%
          add_model(mod)

        mixture_param <-
          parameters(penalty(trans = NULL), mixture())

        glmn_grid <- grid_regular(mixture_param, levels = c(10, 10))


        ctrl <- control_grid(save_pred = TRUE, verbose = F)

        glmn_tune <-
          tune_grid(
            wfl,
            resamples = cv_splits,
            grid = glmn_grid,
            metrics = metric_set(rmse),
            control = ctrl
          )


        best_parameters <- select_best(glmn_tune, metric = "rmse")

        wf1_final <-
          wfl %>%
          finalize_workflow(best_parameters)


        get_model <- function(x) {
          pull_workflow_fit(x) %>% tidy()
        }
        ctrl <- control_resamples(extract = get_model)

        glmn_cv_final <-
          tune::fit_resamples(wfl_final, cv_splits, control = ctrl)

        all_coef <- map_dfr(glmn_cv_final$.extracts, ~ .x[[1]][[1]])

        all_coef_avg <- all_coef %>% group_by(term) %>% summarise(cv_mean = mean(estimate))

        all_coef_avg$abs_cv_mean <- abs(all_coef_avg$cv_mean)
        all_coef_avg$environment <- environment


        return(all_coef_avg)

      }


      all_envs = METData$info_environments$IDenv


      res_all_envs = lapply(
        all_envs,
        FUN = function(x) {
          marker_effect_per_env(
            geno = METData$geno,
            pheno = METData$pheno,
            trait = trait,
            environment = x
          )
        }
      )

      ## Create a specification of the model before fitting





      METData <- list(
        'geno' = geno,
        'map_markers' = map,
        'pheno' = pheno,
        'compute_ECs' = compute_ECs,
        'env_data' = env_data,
        'info_environments' = info_environments,
        'unique_EC_by_geno' = unique_EC_by_geno,
        'filtering_markers' = filtering_markers,
        'filtering_done' = TRUE
      )




}
