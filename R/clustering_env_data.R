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
    
    options(ggrepel.max.overlaps = Inf)
    set.seed(6)
    
    ## Plot based on weather ECs solely
    
    if (!is.null(path_plots)) {
      path_plots <- file.path(path_plots, 'clustering_analysis')
      if (!dir.exists(path_plots)) {
        dir.create(path_plots, recursive = T)
      }
    }
    
    if (!is.null(weather_ECs)) {
      row.names(weather_ECs) <- weather_ECs$IDenv
      
      weather_ECs_unique <-
        weather_ECs %>% dplyr::select(-IDenv,-year,-location)
      
      cols <-
        names(which(apply(weather_ECs_unique, 2, var) != 0))
      
      weather_ECs_unique <-
        as.data.frame(unique(weather_ECs_unique[, cols]))
      
      
      k <- c(2:(nrow(weather_ECs_unique) - 1))
      if (max(k) > 10) {
        k <- c(2:10)
      }
      
      # Directory for clustering analyses based on weather data only.
      path_plots_w <- file.path(path_plots, 'only_weather')
      if (!dir.exists(path_plots_w)) {
        dir.create(path_plots_w, recursive = T)
      }
      
      evaluate_k_kmeans <- function(k) {
        kclust <- kmeans(weather_ECs_unique,
                         centers = k,
                         nstart = 25)
        
        # First metric: Elbow method: get the percentage of variance explained as a function of the number of clusters
        # Score is the total within-clusters sum of squares
        ss_score <- kclust$tot.withinss
        
        # Second metric: Silhouette score
        sil <-
          cluster::silhouette(kclust$cluster, dist(weather_ECs_unique))
        sil_score <- mean(sil[, 3])
        
        # Plots output for the range of k values
        ## OUTPUT plots: see how environments cluster and which weather-based
        ## covariates might drive the clustering procedure based on PCA
        
        factoextra::fviz_cluster(kclust, data = weather_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots_w,
            paste0(
              'climate_variables_only_clusters_environments_',
              k,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <-
          FactoMineR::PCA(weather_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(
            path_plots_w,
            paste0('PCA_climate_variables_', k, '.pdf')
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        
        return(list("ss_score" = ss_score,
                    "sil_score" = sil_score))
        
      }
      
      metrics_scores <- lapply(k, evaluate_k_kmeans)
      names(metrics_scores) <- k
      
      pdf(file.path(path_plots_w,
                    "plot_sil_score.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['sil_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Average Silhouette Scores',
        frame = FALSE
      )
      dev.off()
      
      pdf(file.path(path_plots_w,
                    "plot_elbow_method.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['ss_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Total within-clusters sum of squares',
        frame = FALSE
      )
      dev.off()
      
      
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
      if (max(k) > 10) {
        k <- c(2:10)
      }
      
      # Directory for clustering analyses based on weather data only.
      path_plots_s <- file.path(path_plots, 'only_soil')
      if (!dir.exists(path_plots_s)) {
        dir.create(path_plots_s, recursive = T)
      }
      
      evaluate_k_kmeans <- function(k) {
        kclust <- kmeans(soil_ECs_unique,
                         centers = k,
                         nstart = 25)
        
        # First metric: Elbow method: get the percentage of variance explained as a function of the number of clusters
        # Score is the total within-clusters sum of squares
        ss_score <- kclust$tot.withinss
        
        # Second metric: Silhouette score
        sil <-
          cluster::silhouette(kclust$cluster, dist(soil_ECs_unique))
        sil_score <- mean(sil[, 3])
        
        # Plots output for the range of k values
        ## OUTPUT plots: see how environments cluster and which weather-based
        ## covariates might drive the clustering procedure based on PCA
        
        factoextra::fviz_cluster(kclust, data = soil_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots_s,
            paste0(
              'climate_variables_only_clusters_environments_',
              k,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <-
          FactoMineR::PCA(soil_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(
            path_plots_s,
            paste0('PCA_climate_variables_', k, '.pdf')
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        
        return(list("ss_score" = ss_score,
                    "sil_score" = sil_score))
        
      }
      
      metrics_scores <- lapply(k, evaluate_k_kmeans)
      names(metrics_scores) <- k
      
      pdf(file.path(path_plots_s,
                    "plot_sil_score.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['sil_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Average Silhouette Scores',
        frame = FALSE
      )
      dev.off()
      
      pdf(file.path(path_plots_s,
                    "plot_elbow_method.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['ss_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Total within-clusters sum of squares',
        frame = FALSE
      )
      dev.off()
      
      
    }
    
    ## Plot based on weather+soil variables together
    
    if (!is.null(soil_ECs) & !is.null(weather_ECs)) {
      all_ECs <-
        merge(soil_ECs,
              weather_ECs %>% dplyr::select(-year,-location),
              by = "IDenv")
      
      row.names(all_ECs) <- all_ECs$IDenv
      
      all_ECs_unique <-
        all_ECs %>% dplyr::select(-IDenv,-year,-location)
      
      cols <-
        names(which(apply(all_ECs_unique, 2, var) != 0))
      
      all_ECs_unique <-
        as.data.frame(unique(all_ECs_unique[, cols]))
      
      k <- c(1:(nrow(all_ECs_unique) - 1))
      if (max(k) > 10) {
        k <- c(2:10)
      }
      
      # Directory for clustering analyses based on weather data only.
      path_plots_all <- file.path(path_plots, 'all_env_data')
      if (!dir.exists(path_plots_all)) {
        dir.create(path_plots_all, recursive = T)
      }
      
      evaluate_k_kmeans <- function(k) {
        kclust <- kmeans(all_ECs_unique,
                         centers = k,
                         nstart = 25)
        
        # First metric: Elbow method: get the percentage of variance explained as a function of the number of clusters
        # Score is the total within-clusters sum of squares
        ss_score <- kclust$tot.withinss
        
        # Second metric: Silhouette score
        sil <-
          cluster::silhouette(kclust$cluster, dist(all_ECs_unique))
        sil_score <- mean(sil[, 3])
        
        # Plots output for the range of k values
        ## OUTPUT plots: see how environments cluster and which weather-based
        ## covariates might drive the clustering procedure based on PCA
        
        factoextra::fviz_cluster(kclust, data = all_ECs_unique, labelsize = 12) +
          theme(axis.text.x = element_text(size = 15),
                title = element_text(size = 15))
        ggsave(
          filename = file.path(
            path_plots_all,
            paste0(
              'climate_variables_only_clusters_environments_',
              k,
              '.pdf'
            )
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        res.pca <-
          FactoMineR::PCA(all_ECs_unique,  graph = FALSE)
        factoextra::fviz_pca_biplot(res.pca, repel = T)
        ggsave(
          filename = file.path(
            path_plots_all,
            paste0('PCA_climate_variables_', k, '.pdf')
          ),
          device = 'pdf',
          height = 8,
          width = 12
        )
        
        return(list("ss_score" = ss_score,
                    "sil_score" = sil_score))
        
      }
      
      metrics_scores <- lapply(k, evaluate_k_kmeans)
      names(metrics_scores) <- k
      
      pdf(file.path(path_plots_all,
                    "plot_sil_score.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['sil_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Average Silhouette Scores',
        frame = FALSE
      )
      dev.off()
      
      pdf(file.path(path_plots_all,
                    "plot_elbow_method.pdf"))
      plot(
        k,
        type = 'b',
        unlist(lapply(metrics_scores, function(x) {
          x[['ss_score']]
        })),
        xlab = 'Number of clusters',
        ylab = 'Total within-clusters sum of squares',
        frame = FALSE
      )
      dev.off()
      
      
    }
  }
