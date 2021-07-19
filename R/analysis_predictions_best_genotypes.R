analysis_predictions_best_genotypes <-
  function(met_pred,
           path_save_results,
           cluster_years = TRUE,
           env_predictors_for_clustering = NULL,
           trait,
           K = 3) {
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
      
      table_res <- dplyr::arrange(table_res,-Mean, sd)
      
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
    ## STEP 2 : categorize environments included in predictions paccording to
    ## weather data and determine best genotypes according to the cluster type.
    ########################################################################
    
    if (cluster_years) {
      if (is.null(env_predictors_for_clustering)) {
        data_to_cluster <-
          met_pred$list_results[[1]]$test %>% dplyr::select(-matches('PC'))
        data_to_cluster <-
          data_to_cluster %>% dplyr::select(-all_of(trait))
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
      
      kclust <- kmeans(data_to_cluster, centers = K)
      
      ## OUTPUT plots: see how environments cluster and which environmental
      ## covariates might drive the clustering procedure based on PCA
      
      factoextra::fviz_cluster(kclust, data = data_to_cluster, labelsize = 7)
      ggsave(
        file.path(path_save_results, 'clusters_env_predicted.pdf'),
        height = 8,
        width = 12
      )
      res.pca <- FactoMineR::PCA(data_to_cluster,  graph = FALSE)
      fviz_pca_biplot(res.pca, repel = T)
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
      
      best_geno_means <- lapply(list_by_cluster, function(x) {
        table_res <- as.data.frame(
          as.data.frame(x) %>%
            group_by(geno_ID) %>%
            dplyr::mutate(Mean = mean(.pred, na.rm = TRUE)) %>%
            dplyr::mutate(sd = sd(.pred, na.rm = TRUE))
        )
        
        table_res <-
          unique(table_res[, c('geno_ID', 'Mean', 'sd')])
        
        table_res <- dplyr::arrange(table_res,-Mean, sd)
        
        return(table_res)
        
      })
      
      names(best_geno_means) <-
        lapply(list_by_cluster, function(x) {
          unique(as.character(x$cluster))
        })
      
      saveRDS(
        object = best_geno_means,
        file = file.path(path_save_results, 'geno_means_by_cluster.RDS')
      )
      
      
    }
    
    
    
    
    
  }
