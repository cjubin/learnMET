#' Clustering of environments solely based on environmental information
#' @description
#' Clustering of environments using a K-means algorithm based on:
#' \enumerate{
#' \item Only climate variables
#' \item Only soil variables
#' \item All environmental data together
#' }
#' @param weather_ECs a \code{data.frame} with climate-based covariates (predictor variables)
#'   with three columns year, location and IDenv.
#'   Default can be `NULL`, if no weather-based data available/should be used.
#' @param soil_ECs a \code{data.frame} with soil-based covariates (predictor variables)
#'   with three columns year, location and IDenv.
#'   Default can be `NULL`, if no soil data available.
#' @param k \code{numeric} Number of clusters to use.
#' @param path_plots \code{character} Path where clustering plots should be saved.
#' @return
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'



clustering_env_data <-
  function(weather_ECs = NULL,
           soil_ECs = NULL,
           path_plots = NULL) {
    ## Plot based on weather ECs solely
    
    if (!is.null(weather_ECs)) {
      row.names(weather_ECs) <- weather_ECs$IDenv
      
      weather_ECs_unique <-
        weather_ECs %>% dplyr::select(-IDenv,-year,-location)
      
      cols <-
        names(which(apply(weather_ECs_unique, 2, var) != 0))
      
      weather_ECs_unique <-
        as.data.frame(unique(weather_ECs_unique[, cols]))
      
      
      k <- c(1:(nrow(weather_ECs_unique) - 1))
      if (max(k) > 8) {
        k <- c(1:8)
      }
      
      for (j in 1:length(k)) {
        K <- k[j]
        kclust <- kmeans(weather_ECs_unique, centers = K)
        
        
        
        
        ## OUTPUT plots: see how environments cluster and which environmental
        ## covariates might drive the clustering procedure based on PCA
        
        
        factoextra::fviz_cluster(kclust, data = weather_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots,
            paste0(
              'climate_variables_only_clusters_environments_',
              K,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <- FactoMineR::PCA(weather_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(path_plots, paste0(
            'PCA_climate_variables_', K, '.pdf'
          )),
          device = 'pdf',
          height = 8,
          width = 12
        )
        
        
      }
    }
    
    
    ## Plot based on soil ECs solely
    
    if (!is.null(soil_ECs)) {
      row.names(soil_ECs) <- soil_ECs$IDenv
      
      soil_ECs_unique <-
        soil_ECs %>% dplyr::select(-IDenv,-year,-location)
      
      cols <-
        names(which(apply(soil_ECs_unique, 2, var) != 0))
      
      soil_ECs_unique <-
        as.data.frame(unique(soil_ECs_unique[, cols]))
      
      k <- c(1:(nrow(soil_ECs_unique) - 1))
      if (max(k) > 8) {
        k <- c(1:8)
      }
      
      for (j in 1:length(k)) {
        K <- k[j]
        kclust <- kmeans(soil_ECs_unique, centers = K)
        
        
        
        
        ## OUTPUT plots: see how environments cluster and which environmental
        ## covariates might drive the clustering procedure based on PCA
        
        
        factoextra::fviz_cluster(kclust, data = soil_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots,
            paste0(
              'soil_variables_only_clusters_environments_',
              K,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <- FactoMineR::PCA(soil_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(path_plots, paste0('PCA_soil_variables_', K, '.pdf')),
          device = 'pdf',
          height = 8,
          width = 12
        )
        
        
      }
      
    }
    
    ## Plot based on weather+soil variables together
    
    if (!is.null(soil_ECs) & !is.null(weather_ECs)) {
      
      all_ECs <- merge(soil_ECs, weather_ECs %>% dplyr::select(-year,-location), by = "IDenv")
      
      row.names(all_ECs) <- all_ECs$IDenv
      
      all_ECs_unique <-
        all_ECs %>% dplyr::select(-IDenv,-year,-location)
      
      cols <-
        names(which(apply(all_ECs_unique, 2, var) != 0))
      
      all_ECs_unique <-
        as.data.frame(unique(all_ECs_unique[, cols]))
      
      k <- c(1:(nrow(all_ECs_unique) - 1))
      if (max(k) > 8) {
        k <- c(1:8)
      }
      
      for (j in 1:length(k)) {
        K <- k[j]
        kclust <- kmeans(all_ECs_unique, centers = K)
        
        
        
        
        ## OUTPUT plots: see how environments cluster and which environmental
        ## covariates might drive the clustering procedure based on PCA
        
        
        factoextra::fviz_cluster(kclust, data = all_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots,
            paste0(
              'all_env_variables_only_clusters_environments_',
              K,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <- FactoMineR::PCA(all_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(path_plots, paste0(
            'PCA_all_env_variables_', K, '.pdf'
          )),
          device = 'pdf',
          height = 8,
          width = 12
        )
      }
      
    }
  }
