#default elasticnet
#but add rqtl + ridge reg ,and RF/xgboost later


select_markers <- function(METData,
                           trait,
                           method = 'elasticnet',
                           folds_cv = 5,
                           reps = 1) {
  ##

  if (METData$filtering_markers == TRUE &&
      dim(METData$geno[2] < 1000) {
        stop('The number of markers is low and does not need to be further reduced.')
      }



      ## Separate analysis per environment

      marker_effect_per_env <- function(environment, geno, trait) {

        pheno <- METData$pheno[METData$pheno$IDenv == environment,]
        geno <- as.data.frame(METData$geno)
        list_predictors <- colnames(geno)
        geno$geno_ID = row.names(geno)
        pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)

        # Reduce pheno data.frame to trait +  markers as columns

        pheno <- pheno[,c(trait,list_predictors)]

        ## Create the cross-validation randomly splits

        cv_splits <-
          rsample::vfold_cv(pheno, folds = folds_cv, repeats = reps)

        mod <- parsnip::linear_reg(penalty = tune(),
                                   mixture = tune()) %>%
          set_engine("glmnet")



        rec <- recipe(pheno[,1]~., data = pheno) %>%
          update_role(trait,'outcome')
          update_role()
          step_nzv(all_numeric()) %>%
          step_normalize(all_numeric())

      }

      all_envs = METData$info_environments$IDenv

      ## Create a specification of the model before fitting





      METData <- list(
        'geno' = geno,
        'map_markers' = map,
        'pheno' = pheno,
        'compute_ECs' = compute_ECs,
        'env_data' = env_data,
        'info_environments' = info_environments,
        'unique_EC_by_geno' = unique_EC_by_geno,
        'filtering_markers' = filtering_markers
        'filtering_done' = TRUE
      )




}
