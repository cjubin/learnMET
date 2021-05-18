



get_splits_processed_with_method <- function(splits,
                                             method_processing,
                                             trait,
                                             geno_data,
                                             geno_information,
                                             use_selected_markers,
                                             SNPs,
                                             list_env_predictors,
                                             include_env_predictors,
                                             lat_lon_included,
                                             year_included) {
  checkmate::assert_class(splits,
                          "cv_object")
  checkmate::assert_choice(method_processing,
                           choices = c("xgb_ordinal", "xgb_reg","svm_stacking_reg"))
  switch_method <- function(split,
                            trait,
                            geno_data,
                            geno_information,
                            use_selected_markers,
                            SNPs,
                            list_env_predictors,
                            include_env_predictors,
                            lat_lon_included,
                            year_included) {
    switch(
      method_processing,
      xgb_ordinal = xgb_ordinal(
        split = split,
        trait = trait,
        geno_data = geno_data,
        geno_information = geno_information,
        use_selected_markers = use_selected_markers,
        SNPs = use_selected_markers,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included
      ),
      xgb_regression = xgb_regression(
        split = split,
        trait = trait,
        geno_data = geno_data,
        geno_information = geno_information,
        use_selected_markers = use_selected_markers,
        SNPs = use_selected_markers,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included
      ),
      svm_stacking_reg = svm_stacking_reg(
        split = split,
        trait = trait,
        geno_data = geno_data,
        geno_information = geno_information,
        use_selected_markers = use_selected_markers,
        SNPs = use_selected_markers,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included
      )
    )
  }
  
  all_processed_splits <-
    lapply(splits, function(x) {
      switch_method(
        split = x,
        trait = trait,
        geno_data = geno_data,
        geno_information = geno_information,
        use_selected_markers = use_selected_markers,
        SNPs = use_selected_markers,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included
      )
    })
  
  return(all_processed_splits)
}