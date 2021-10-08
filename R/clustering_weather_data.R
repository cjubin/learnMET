#' Compute ECs based on a fixed number of day-windows (fixed number across
#' all environments).
#'
#' @description
#' Compute the environmental covariates based on the daily weather
#' table of an environment (Year x Location), and over a fixed number of time
#' windows, which is common across all environments. The length of day-windows
#' in days in each environment is determined by dividing the total length of
#' the growing season of this environment by the number of windows to use.
#' Each EC is then computed over a fixed number of days within each environment,
#' but the length of the windows can vary across environments.
#'
#' @param table_daily_W \code{data.frame} Object returned by the function
#'   [get_daily_tables_per_env()]
#'
#' @param nb_windows_intervals \code{numeric} Number of day-windows covering
#'   the growing season length (common number of day-windows across all
#'   environments).
#'
#' @param base_temperature \code{numeric} Base temperature (crop growth assumed
#'   to be null below this value.). Default is 10.
#'
#' @param method_GDD_calculation \code{character} Method used to compute the
#'   GDD value, with one out of \code{method_a} or \code{method_b}. \cr
#'   \code{method_a}: No change of the value of \eqn{T_{min}}.
#'   GDD = \eqn{max (\frac{T_{min}+T_{max}}{2} - T_{base},0)}. \cr
#'   \code{method_b}: If \eqn{T_{min}} < \eqn{T_{base}}, change \eqn{T_{min}}
#'   to \eqn{T_{min}} = \eqn{T_{base}}. \cr
#'   Default = \code{method_b}.
#'
#' @return An object of class \code{data.frame} with
#'   10 x number_total_fixed_windows + 1 last column (IDenv):
#'   \enumerate{
#'     \item mean_TMIN: number_total_fixed_windows columns, indicating the
#'     average minimal temperature over the respective day-window.
#'     \item mean_TMAX: number_total_fixed_windows columns, indicating the
#'     average maximal temperature over the respective day-window.
#'     \item mean_TMEAN: number_total_fixed_windows columns, indicating the
#'     average mean temperature over the respective day-window.
#'     \item freq_TMAX_sup30: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 30°C over the respective
#'     day-window.
#'     \item freq_TMAX_sup35: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 35°C over the respective
#'     day-window.
#'     \item sum_GDD: number_total_fixed_windows columns, indicating the
#'     growing degree days over the respective day-window.
#'     \item sum_PTT: number_total_fixed_windows columns, indicating the
#'     accumulated photothermal time over the respective day-window.
#'     \item sum_P: number_total_fixed_windows columns, indicating the
#'     accumulated precipitation over the respective day-window.
#'     \item freq_P_sup10: number_total_fixed_windows columns, indicating the
#'     frequency of days with total precipitation superior to 10 mm over the
#'     respective day-window.
#'     \item sum_solar_radiation: number_total_fixed_windows columns, indicating
#'     the accumulated incoming solar radiation over the respective day-window.
#'     \item IDenv \code{character} ID of the environment (Location_Year)
#'    }
#'   @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#'   @export
#'



clustering_weather_data <- function(weather_ECs, k = c(1:7),path_plots) {
  weather_ECs <-
    weather_ECs %>% dplyr::select(-IDenv, -year, -location)
  
  cols <-
    names(which(apply(weather_ECs, 2, var) != 0))
  
  weather_ECs <-
    as.data.frame(unique(weather_ECs[, cols]))
  
  for (j in 1:length(k)) {
    K <- k[j]
    kclust <- kmeans(weather_ECs, centers = K)
    
    
    
    
    ## OUTPUT plots: see how environments cluster and which environmental
    ## covariates might drive the clustering procedure based on PCA
    
    
    factoextra::fviz_cluster(kclust, data = weather_ECs, labelsize = 12) +theme(axis.text.x = element_text(size = 15), title = element_text(size = 15))
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
