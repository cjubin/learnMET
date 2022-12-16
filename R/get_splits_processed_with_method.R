#' Attribute a processing method for each list of training/test splits
#'
#' @description
#' Attribute the function for processing the dataset according to the processing parameters.
#'
#' @param splits an object of class `cv_object`
#' @param prediction_method
#'
#' @return a list of objects of the class of the chosen prediction method
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @export
#'
#'



get_splits_processed_with_method <- function(splits,
                                             prediction_method,
                                             trait,
                                             geno,
                                             env_predictors,
                                             info_environments,
                                             use_selected_markers,
                                             SNPs,
                                             list_env_predictors,
                                             include_env_predictors,
                                             lat_lon_included,
                                             year_included,
                                             ...) {
  switch_method <- function(split,
                            prediction_method,
                            trait,
                            geno,
                            env_predictors,
                            info_environments,
                            use_selected_markers,
                            SNPs,
                            list_env_predictors,
                            include_env_predictors,
                            lat_lon_included,
                            year_included,
                            ...) {
    switch(
      prediction_method,
      xgb_reg_1 = xgb_reg_1(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      xgb_reg_2 = xgb_reg_2(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      xgb_reg_3 = xgb_reg_3(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      rf_reg_1 = rf_reg_1(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      rf_reg_2 = rf_reg_2(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      rf_reg_3 = rf_reg_3(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      DL_reg_1 = DL_reg_1(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      DL_reg_2 = DL_reg_2(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      DL_reg_3 = DL_reg_3(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      stacking_reg_1 = stacking_reg_1(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      stacking_reg_2 = stacking_reg_2(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      ),
      stacking_reg_3 = stacking_reg_3(
        split = split,
        trait = trait,
        geno = geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        use_selected_markers = use_selected_markers,
        SNPs = SNPs,
        list_env_predictors = list_env_predictors,
        include_env_predictors = include_env_predictors,
        lat_lon_included = lat_lon_included,
        year_included = year_included,
        ...
      )
    )
  }
  
  
  all_processed_splits = list()
  length(all_processed_splits) <- length(splits)
  
  optional_args <- as.list(match.call(expand.dots = TRUE))
  optional_args <- optional_args[optional_args != "splits"]
  optional_args <-
    optional_args[optional_args != "get_splits_processed_with_method"]
  for (i in 1:length(all_processed_splits)) {
    optional_args$split <- splits[[i]]
    all_processed_splits[[i]] <-
      do.call(switch_method, args = optional_args)
  }
  
  
  return(all_processed_splits)
}