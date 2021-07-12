#' Attribute 
#'
#' @description
#' Implement trait prediction based on SNP and environmental data
#' with selection of prediction methods among Machine Learning approaches.
#'
#' This function should be used to assess the predictive ability according to
#' a cross-validation scheme determined by the user.
#'
#' @param splits
#' @param method_processing
#'
#' @return 
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' 
#' 



get_splits_processed_with_method <- function(splits,
                                             method_processing,
                                             trait,
                                             geno_data,
                                             env_predictors,
                                             info_environments,
                                             unique_EC_by_geno,
                                             geno_information,
                                             use_selected_markers,
                                             SNPs,
                                             list_env_predictors,
                                             include_env_predictors,
                                             lat_lon_included,
                                             year_included) {
  
  switch_method <- function(split,
                            method_processing,
                            trait,
                            geno_data,
                            env_predictors,
                            info_environments,
                            unique_EC_by_geno,
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
        split=split,
        trait=trait,
        geno_data=geno_data,
        env_predictors = env_predictors,
        info_environments = info_environments,
        unique_EC_by_geno=unique_EC_by_geno,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included
      ),
      xgb_reg = xgb_reg(
        split=split,
        trait=trait,
        geno_data=geno_data,
        env_predictors = env_predictors,
        info_environments = info_environments,
        unique_EC_by_geno=unique_EC_by_geno,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included
      ),
      DL_reg = DL_reg(
        split=split,
        trait=trait,
        geno_data=geno_data,
        env_predictors = env_predictors,
        info_environments = info_environments,
        unique_EC_by_geno=unique_EC_by_geno,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included
      ),
      svm_stacking_reg = svm_stacking_reg(
        split=split,
        trait=trait,
        geno_data=geno_data,
        env_predictors = env_predictors,
        info_environments = info_environments,
        unique_EC_by_geno=unique_EC_by_geno,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included
      )
    )
  }
  
  all_processed_splits <-
    lapply(splits, function(x) {
      switch_method(
        split = x,
        method_processing = method_processing, 
        trait=trait,
        geno_data=geno_data,
        env_predictors = env_predictors,
        info_environments = info_environments,
        unique_EC_by_geno=unique_EC_by_geno,
        geno_information=geno_information,
        use_selected_markers=use_selected_markers,
        SNPs=SNPs,
        list_env_predictors=list_env_predictors,
        include_env_predictors=include_env_predictors,
        lat_lon_included=lat_lon_included,
        year_included=year_included
      )
    })
  
  return(all_processed_splits)
}