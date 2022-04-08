#' Analysis of prediction results from [predict_trait_MET()] by location and by
#' environmental cluster.
#'
#' @description
#' This function analyses at each individual location the best genotypes based
#' on their predicted performance over years. List of mean genotype
#' performance (simple average based on predictions) at each location are saved.
#' Second part of the function performs a K-means clustering of environments
#' based on the env. covariables characterizing each of these environments.
#' List of mean genotype performance (based on predictions) for each
#' environmental cluster are saved.
#' No output from this function, all results are direcly saved based on the
#' `path_save_results` argument.
#'
#' @param met_pred Output from [predict_trait_MET()] function of class
#'   `met_pred`.
#'
#' @param path_save_results Path where plots and results should be saved.
#'
#' @param cluster_envs \code{logical} indicates if k-means algorithm
#'   on environmental data should be done. Default is `TRUE`.
#'
#' @param env_predictors_for_clustering \code{character} Vector giving the
#'   names of environmental covariates to be used in the K-means clustering
#'   analysis.
#'
#' @param K \code{integer} number of clusters to use in K-means algorithm
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export


analysis_predictions <-
  function(met_pred,
           path_save_results,
           cluster_envs = TRUE,
           env_predictors_for_clustering = NULL,
           name_trait = NULL,
           K = 4) {
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
          group_by(geno_ID) %>%
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
    
    
    ########################################################################
    ## STEP 2 : categorize environments included in predictions according to
    ## weather data and determine best genotypes according to the cluster type.
    ########################################################################
    
    
    if (cluster_envs) {
      if (is.null(env_predictors_for_clustering)) {
        # Identify names of weather-based variables
        
        data_to_cluster <- met_pred$list_results[[1]]$test
        
        columns_weather <-
          data_to_cluster %>% dplyr::select(-yld_bu_ac, -year, -location) %>% group_by(IDenv) %>% group_map( ~
                                                                                                               apply(.x, 2, var))
        columns_weather <-
          names(columns_weather[[1]][which(columns_weather[[1]] == 0)])
        
        
        data_to_cluster <-
          met_pred$list_results[[1]]$test %>% dplyr::select(all_of(columns_weather), IDenv)
        
      } else{
        data_to_cluster <-
          met_pred$list_results[[1]]$test %>% dplyr::select(contains(env_predictors_for_clustering))
      }
      
      
      cols <-
        c(names(which(
          apply(data_to_cluster %>% select(-IDenv), 2, var) != 0
        )), 'IDenv')
      
      data_to_cluster <-
        as.data.frame(unique(data_to_cluster[, cols]))
      rownames(data_to_cluster) <- data_to_cluster$IDenv
      data_to_cluster <- data_to_cluster %>% select(-IDenv)
      
      kclust <- kmeans(data_to_cluster, centers = K, nstart = 25)
      
      ## OUTPUT plots: see how environments cluster and which environmental
      ## covariates might drive the clustering procedure based on PCA
      
      factoextra::fviz_cluster(kclust, data = data_to_cluster, labelsize = 7)
      ggsave(
        file.path(path_save_results, 'clusters_env_predicted.pdf'),
        height = 8,
        width = 12
      )
      res.pca <- FactoMineR::PCA(data_to_cluster,  graph = FALSE)
      factoextra::fviz_pca_biplot(res.pca, repel = T)
      ggsave(
        file.path(path_save_results, 'PCA_env_predicted.pdf'),
        height = 8,
        width = 12
      )
      
      res_kclust <- data.frame(kclust$cluster)
      res_kclust$IDenv <- row.names(res_kclust)
      
      ## Create tables of mean performance of genotypes within each cluster
      
      pred_data <- met_pred$list_results[[1]]$predictions_df
      pred_data$cluster <-
        res_kclust[match(pred_data$IDenv, res_kclust$IDenv), 'kclust.cluster']
      
      list_by_cluster <- split(pred_data, pred_data$cluster)
      
      best_geno_means_cluster <-
        lapply(list_by_cluster, function(x) {
          table_res <- as.data.frame(
            as.data.frame(x) %>%
              group_by(geno_ID) %>%
              dplyr::mutate(Mean = mean(.pred, na.rm = TRUE)) %>%
              dplyr::mutate(sd = sd(.pred, na.rm = TRUE))
          )
          
          table_res <-
            unique(table_res[, c('geno_ID', 'Mean', 'sd')])
          
          table_res <- dplyr::arrange(table_res, -Mean, sd)
          
          return(table_res)
          
        })
      
      names(best_geno_means_cluster) <-
        lapply(list_by_cluster, function(x) {
          unique(as.character(x$cluster))
        })
      
      saveRDS(
        object = best_geno_means_cluster,
        file = file.path(path_save_results, 'geno_means_by_cluster.RDS')
      )
      
      
    }
    
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
    
    if (is.null(name_trait)){
      name_trait <- 'trait'
    }
      
    
    ggplot(top_by_loc_pred, aes(x = geno_ID, y = .pred)) + geom_boxplot() + facet_wrap( ~ location, scales = 'free_x') +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))+
      labs(title = paste0('Predicted distribution of ',trait,' for each genotype–location combination across ',length(unique(pred_data$year)), ' years'))
    
    ggsave(file.path(
      path_save_results,
      'top_predicted_by_loc_across_years.pdf'
    ))
    
    
    
    
  }
