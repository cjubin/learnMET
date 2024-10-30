#' Compute the DTW distance using the dtwcluster R package & compute features
#'
#' @description
#' Compute DTW distance using the dtw_basic function from dtwcluster R package
#' and  multi-dimensional scaling on the obtained distance to get features for
#' ML analyses
#'
#' @param a list of data.frame with the following columns extracted
#' from NASA POWER data, according to requested parameters:
#' \enumerate{
#'   \item longitude \code{numeric}
#'   \item latitude \code{numeric}
#'   \item YEAR \code{numeric}
#'   \item MM \code{integer}
#'   \item DD \code{integer}
#'   \item DOY \code{integer}
#'   \item YYYYMMDD \code{Date}
#'   \item RH2M \code{numeric}
#'   \item T2M \code{numeric}
#'   \item T2M_MIN \code{numeric}
#'   \item T2M_MAX \code{numeric}
#'   \item PRECTOTCORR \code{numeric}
#'   \item ALLSKY_SFC_SW_DWN \code{numeric}
#'   \item T2MDEW \code{numeric}
#'   \item IDenv \code{character} ID environment for which weather data were
#'   downloaded.
#'   \item length.gs \code{difftime} length in days of the growing season
#'   for the environment.
#'   }
#' 
#' @return a data.frame \code{data.frame} with the results from the multi
#' dimensional scaling analysis.
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#' @export



dtw_distance <-
  function(daily_weather_data,
           path_data,
           reduction_dimensionality = FALSE,
           ...) {
      if (!is.null(path_data)) {
          path_data <- file.path(path_data, 'plots_dtw')
          if (!dir.exists(path_data)) {
              dir.create(path_data, recursive = T)
          }
      } else{
      stop("Please give the path to save the plots.\n")
      }
      
      cat("Calculation of dtw distance\n")
      for (l in 1:length(daily_weather_data)){
          daily_weather_data[[l]] <- as.matrix(as.data.frame(daily_weather_data[[l]]  %>%
                                                             dplyr::select(RH2M,T2M,T2M_MIN,T2M_MAX,PRECTOTCORR)))
      }
      distance_dtw <- as.matrix(proxy::dist(
          daily_weather_data,
          method = "dtw_basic",
          norm = "L1",
          step.pattern=symmetric1,
          normalize = F
        ))
      
      if (reduction_dimensionality){
          res_cmdscale <- as.data.frame(cmdscale(distance_dtw))
          res_cmdscale$IDenv <- row.names(res_cmdscale)
          pcoa <- cmdscale(distance_dtw)
      } else{
          distance_dtw <- as.data.frame(distance_dtw)
          return(distance_dtw)
          
      }
      
      #ordiplot(pcoa, display = 'sites', type = 'text')
      #ggsave(
      #    filename = file.path(
      #      path_data,
      #      paste0(
      #        'dtw_plot.png'
      #      )
      #    ),
      #    device = 'png',
      #    height = 8,
      #    width = 12
       # )
        #d_norm=(distance_dtw-min(distance_dtw))/(max(distance_dtw)-min(distance_dtw))
        #sim=1-d_norm
        #colnames(sim)<-row.names(sim)<-colnames(sim)

      return(res_cmdscale)

}