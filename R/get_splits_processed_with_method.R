#' Attribute 
#'
#' @description
#' Implement trait prediction based on SNP and environmental data
#' with selection of prediction methods among Machine Learning approaches.
#'
#' This function should be used to assess the predictive ability according to
#' a cross-validation scheme determined by the user.
#'
#' @param splits zw4zw
#' @param prediction_method rze
#'
#' @return jj
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' 
#' 



get_splits_processed_with_method <- function(splits,
                                             prediction_method,
                                             trait,
                                             geno,
                                             env_predictors,
                                             info_environments,
                                             geno_information,
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
                            geno_information,
                            use_selected_markers,
                            SNPs,
                            list_env_predictors,
                            include_env_predictors,
                            lat_lon_included,
                            year_included,
                            ...) {
    switch(
      prediction_method,
      xgb_reg = xgb_reg(
        split=split,
        trait=trait,
        geno=geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included,
        ...
      ),
      DL_reg = DL_reg(
        split=split,
        trait=trait,
        geno=geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included,
        ...
      ),
      stacking_reg_1 = stacking_reg_1(
        split=split,
        trait=trait,
        geno=geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included,
        ...
      ),
      stacking_reg_2 = stacking_reg_2(
        split=split,
        trait=trait,
        geno=geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included,
        ...
      ),
      stacking_reg_3 = stacking_reg_3(
        split=split,
        trait=trait,
        geno=geno,
        env_predictors = env_predictors,
        info_environments = info_environments,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included,
        ...
      )
    )
  }
  
  
  all_processed_splits = list()
  length(all_processed_splits) <- length(splits)

  optional_args <- as.list(match.call(expand.dots=TRUE))
  optional_args <- optional_args[optional_args != "splits"]
  
  
  for (i in 1:length(all_processed_splits)) {
    optional_args$split <- splits[[i]]
    all_processed_splits[[i]] <-
      do.call(switch_method, args = optional_args)
  }
  
  
  return(all_processed_splits)
}