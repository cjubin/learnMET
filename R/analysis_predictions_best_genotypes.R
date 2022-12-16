#' Analysis of prediction results from [predict_trait_MET()] by location.
#'
#' @description
#' This function analyses at each individual location the best genotypes based
#' on their predicted performance over years. List of mean genotype
#' performance (simple average based on predictions) at each location are saved.
#'
#' No output from this function, all results are direcly saved based on the
#' `path_save_results` argument.
#'
#' @param met_pred Output from [predict_trait_MET()] function of class
#'   `met_pred`.
#'
#' @param path_save_results Path where plots and results should be saved.
#'
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @export


analysis_predictions <-
  function(met_pred,
           path_save_results,
           trait = NULL) {
    checkmate::assert_class(met_pred, 'met_pred')
    
    # Check the path_folder: create if does not exist
    if (!dir.exists(path_save_results)) {
      dir.create(path_save_results, recursive = T)
    }
    
    
    pred_data <- met_pred$list_results[[1]]$predictions_df
    
    checkmate::test_numeric(pred_data$.pred, unique = FALSE)
    
    ## STEP 1 : summary tables of genotypes based on highest mean and sd, for
    ## each location, across predicted years
    
    list_by_loc <- split(pred_data, pred_data$location)
    
    best_geno_means <- lapply(list_by_loc, function(x) {
      table_res <- as.data.frame(
        as.data.frame(x) %>%
          dplyr::group_by(geno_ID) %>%
          dplyr::mutate(Mean = mean(.pred, na.rm = TRUE)) %>%
          dplyr::mutate(sd = sd(.pred, na.rm = TRUE))
      )
      
      table_res <-
        unique(table_res[, c('geno_ID', 'Mean', 'sd')])
      
      table_res <- dplyr::arrange(table_res, -Mean, sd)
      
      return(table_res)
      
    })
    
    names(best_geno_means) <-
      lapply(list_by_loc, function(x) {
        unique(as.character(x$location))
      })
    
    saveRDS(
      object = best_geno_means,
      file = file.path(path_save_results, 'geno_means_by_location.RDS')
    )
    
    
    
    
    ## Plot by location of top varieties at each location across years
    
    top_by_loc <-
      reshape2::melt(best_geno_means, id.vars = c("Mean", "geno_ID", "sd")) %>% dplyr::group_by(L1)  %>% dplyr::top_n(., n =
                                                                                                                        40, wt = Mean)
    colnames(top_by_loc)[4] <- 'location'
    
    top_by_loc_pred <-
      plyr::join(pred_data,
                 top_by_loc,
                 type = 'left',
                 by = c('geno_ID', 'location'))
    top_by_loc_pred <-
      top_by_loc_pred[!is.na(top_by_loc_pred$Mean), c('.pred', 'year', 'location', 'geno_ID')]
    
    if (is.null(trait)){
      trait <- 'trait'
    }
      
    
    ggplot(top_by_loc_pred, aes(x = geno_ID, y = .pred)) + geom_boxplot() + facet_wrap( ~ location, scales = 'free_x') +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))+
      labs(title = paste0('Predicted distribution of ',trait,' for each genotype by location combination across ',length(unique(pred_data$year)), ' years'))
    
    ggsave(file.path(
      path_save_results,
      'top_predicted_by_loc_across_years.pdf'
    ),width=35,height=20,units = 'cm')
    
    
    
    ggplot(top_by_loc, aes(y = geno_ID, x = location, fill= .pred)) +
      geom_tile() +
      viridis::scale_fill_viridis()
    
    ggsave(file.path(
      path_save_results,
      'heatmap_best_variety_by_location.pdf'
    ),width=30,height=30,units = 'cm')
    
    
    
    
  }
