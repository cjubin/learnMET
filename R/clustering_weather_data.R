clustering_weather_data <- function(weather_ECs, k = c(1:7),path_plots) {
  weather_ECs <-
    weather_ECs %>% select(-IDenv, -year, -location)
  
  cols <-
    names(which(apply(weather_ECs, 2, var) != 0))
  
  weather_ECs <-
    as.data.frame(unique(weather_ECs[, cols]))
  
  for (j in 1:length(k)) {
    K <- k[j]
    kclust <- kmeans(weather_ECs, centers = K)
    
    
    
    
    ## OUTPUT plots: see how environments cluster and which environmental
    ## covariates might drive the clustering procedure based on PCA
    
    
    factoextra::fviz_cluster(kclust, data = weather_ECs, labelsize = 14) +theme(axis.text.x = element_text(size = 15), title = element_text(size = 15))
    ggsave(
      file.path(path_plots, paste0('clusters_env_predicted_',K,'.pdf')),
      height = 8,
      width = 12
    )
    res.pca <- FactoMineR::PCA(weather_ECs,  graph = FALSE)
    factoextra::fviz_pca_biplot(res.pca, repel = T)
    ggsave(
      file.path(path_plots, paste0('PCA_env_predicted_',K,'.pdf')),
      height = 8,
      width = 12
    )
    
    
  }
  
  
}
