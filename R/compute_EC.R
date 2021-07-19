#' Compute ECs based on time windows of fixed length.
#'
#' @description
#' Compute the environmental covariates based on the daily weather
#' table of an environment (Year x Location), and over time windows of fixed
#' length. Each EC is computed over a fixed certain number of days, given by
#' the parameter "duration_time_window_days". The maximum number of time
#' windows (e.g. the total number of ECs) is determined by the parameter
#' number_total_fixed_windows.
#'
#' @param table_daily_W \code{data.frame} Object returned by the function
#'   get_daily_tables_per_env()
#'
#' @param duration_time_window_days \code{numeric} Number of days spanned by a
#'   time window
#'
#' @param length_minimum_gs \code{numeric} Length of the shortest growing season
#'   length. Used to calculate the maximum number of time windows to use
#'   (is determined based on the shortest growing season length).
#'
#' @param base_temperature \code{numeric} Base temperature (crop growth assumed
#'   to be null below this value.) Default is 10.
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
#'     average minimal temperature over the respective time window.
#'     \item mean_TMAX: number_total_fixed_windows columns, indicating the
#'     average maximal temperature over the respective time window.
#'     \item mean_TMEAN: number_total_fixed_windows columns, indicating the
#'     average mean temperature over the respective time window.
#'     \item freq_TMAX_sup30: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 30°C over the respective
#'     time window.
#'     \item freq_TMAX_sup35: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 35°C over the respective
#'     time window.
#'     \item sum_GDD: number_total_fixed_windows columns, indicating the
#'     growing degree days over the respective time window.
#'     \item sum_PTT: number_total_fixed_windows columns, indicating the
#'     accumulated photothermal time over the respective time window.
#'     \item sum_P: number_total_fixed_windows columns, indicating the
#'     accumulated precipitation over the respective time window.
#'     \item freq_P_sup10: number_total_fixed_windows columns, indicating the
#'     frequency of days with total precipitation superior to 10 mm over the
#'     respective time window.
#'     \item sum_solar_radiation: number_total_fixed_windows columns, indicating
#'     the accumulated incoming solar radiation over the respective time window.
#'     \item IDenv \code{character} ID of the environment (Location_Year)
#'    }
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'


compute_EC_fixed_length_window <- function(table_daily_W,
                                           length_minimum_gs,
                                           base_temperature = 10,
                                           method_GDD_calculation =
                                             c('method_b'),
                                           duration_time_window_days = 10,
                                           ...) {
 
  number_total_fixed_windows <-
    floor(length_minimum_gs / duration_time_window_days)
  
  # Calculation GDD
  table_daily_W$TMIN_GDD = table_daily_W$T2M_MIN
  table_daily_W$TMAX_GDD = table_daily_W$T2M_MAX
  
  if (method_GDD_calculation == 'method_b') {
    # Method b: when the minimum temperature T_min is below the T_base:
    # Any temperature below T_base is set to T_base before calculating the
    # average. https://en.wikipedia.org/wiki/Growing_degree-day
    
    table_daily_W$TMIN_GDD[table_daily_W$TMIN_GDD < base_temperature] <-
      base_temperature
    table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD < base_temperature] <-
      base_temperature
  }
  
  # The maximum temperature is usually capped at 30 °C for GDD calculation.
  
  table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD > 30] <- 30
  table_daily_W$TMEAN_GDD <-
    (table_daily_W$TMAX_GDD + table_daily_W$TMIN_GDD) / 2
  table_daily_W$GDD = table_daily_W$TMEAN_GDD - 10
  
  if (method_GDD_calculation == 'method_a') {
    table_daily_W$GDD[table_daily_W$GDD < 0] <- 0
  }
  
  
  # Calculation day length
  
  table_daily_W$day_length <-
    daylength(lat = table_daily_W$LAT, day_of_year = table_daily_W$DOY)
  table_daily_W$PhotothermalTime <-
    table_daily_W$day_length * table_daily_W$GDD
  
  
  mean_TMIN <- zoo::rollapply(table_daily_W$T2M_MIN,
                              width = duration_time_window_days,
                              mean,
                              by = duration_time_window_days)
  
  mean_TMAX = zoo::rollapply(table_daily_W$T2M_MAX,
                             width = duration_time_window_days,
                             mean,
                             by = duration_time_window_days)
  
  
  mean_TMEAN = zoo::rollapply(table_daily_W$T2M,
                              width = duration_time_window_days,
                              mean,
                              by = duration_time_window_days)
  
  freq_TMAX_sup30 = zoo::rollapply(table_daily_W$T2M_MAX,
                                   width = duration_time_window_days,
                                   function(x) {
                                     length(which(x > 30)) / length(x)
                                   },
                                   by = duration_time_window_days)
  
  
  freq_TMAX_sup35 = zoo::rollapply(table_daily_W$T2M_MAX,
                                   width = duration_time_window_days,
                                   function(x) {
                                     length(which(x > 35)) / length(x)
                                   },
                                   by = duration_time_window_days)
  
  
  
  sum_GDD = zoo::rollapply(table_daily_W$GDD,
                           width = duration_time_window_days,
                           sum,
                           by = duration_time_window_days)
  
  sum_PTT = zoo::rollapply(table_daily_W$PhotothermalTime,
                           width = duration_time_window_days,
                           sum,
                           by = duration_time_window_days)
  
  sum_P = zoo::rollapply(table_daily_W$PRECTOT,
                         width = duration_time_window_days,
                         sum,
                         by = duration_time_window_days)
  
  freq_P_sup10 = zoo::rollapply(table_daily_W$PRECTOT,
                                width = duration_time_window_days,
                                function(x) {
                                  length(which(x > 10)) / length(x)
                                },
                                by = duration_time_window_days)
  
  sum_solar_radiation = zoo::rollapply(table_daily_W$ALLSKY_SFC_SW_DWN,
                                       width = duration_time_window_days,
                                       sum,
                                       by = duration_time_window_days)
  
  table_EC <-
    data.frame(
      mean_TMIN,
      mean_TMAX,
      mean_TMEAN,
      freq_TMAX_sup30,
      freq_TMAX_sup35,
      sum_GDD,
      sum_PTT,
      sum_P,
      freq_P_sup10,
      sum_solar_radiation
    )
  
  if (nrow(table_EC) > number_total_fixed_windows) {
    table_EC <- table_EC[1:number_total_fixed_windows, ]
  }
  
  # Format for final EC table per environment
  # Each cell represents the value of the EC for this time window, e.g.
  # represents an EC on its own. Therefore, each cell should represent one
  # column.
  
  grid_tab <-
    as.data.frame(expand.grid(colnames(table_EC), row.names(table_EC)))
  grid_tab <- grid_tab[order(grid_tab$Var1), ]
  row.names(grid_tab) <- NULL
  
  table_EC_long <-
    data.frame(
      t(table_EC$mean_TMIN),
      t(table_EC$mean_TMAX),
      t(table_EC$mean_TMEAN),
      t(table_EC$freq_TMAX_sup30),
      t(table_EC$freq_TMAX_sup35),
      t(table_EC$sum_GDD),
      t(table_EC$sum_PTT),
      t(table_EC$sum_P),
      t(table_EC$freq_P_sup10),
      t(table_EC$sum_solar_radiation)
    )
  
  colnames(table_EC_long) <-
    paste0(grid_tab$Var1, '_', grid_tab$Var2)
  table_EC_long$IDenv <- unique(table_daily_W$IDenv)
  
  return(table_EC_long)
  
}

#' Compute ECs based on a fixed number of time windows.
#'
#' @description
#' Compute the environmental covariates based on the daily weather
#' table of an environment (Year x Location), and over a fixed number of time
#' windows, which is common across all environments. The length of time windows
#' in days in each environment is determined by dividing the total length of
#' the growing season of this environment by the number of windows to use.
#' Each EC is then computed over a fixed number of days, given by
#' the calculated duration_time_window_days.
#'
#' @param table_daily_W \code{data.frame} Object returned by the function
#'   get_daily_tables_per_env()
#'
#' @param nb_windows_intervals \code{numeric} Number of time windows covering
#'   the growing season length (common number of time windows across all
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
#'     average minimal temperature over the respective time window.
#'     \item mean_TMAX: number_total_fixed_windows columns, indicating the
#'     average maximal temperature over the respective time window.
#'     \item mean_TMEAN: number_total_fixed_windows columns, indicating the
#'     average mean temperature over the respective time window.
#'     \item freq_TMAX_sup30: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 30°C over the respective
#'     time window.
#'     \item freq_TMAX_sup35: number_total_fixed_windows columns, indicating the
#'     frequency of days with maximum temperature over 35°C over the respective
#'     time window.
#'     \item sum_GDD: number_total_fixed_windows columns, indicating the
#'     growing degree days over the respective time window.
#'     \item sum_PTT: number_total_fixed_windows columns, indicating the
#'     accumulated photothermal time over the respective time window.
#'     \item sum_P: number_total_fixed_windows columns, indicating the
#'     accumulated precipitation over the respective time window.
#'     \item freq_P_sup10: number_total_fixed_windows columns, indicating the
#'     frequency of days with total precipitation superior to 10 mm over the
#'     respective time window.
#'     \item sum_solar_radiation: number_total_fixed_windows columns, indicating
#'     the accumulated incoming solar radiation over the respective time window.
#'     \item IDenv \code{character} ID of the environment (Location_Year)
#'    }
#'   @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#'   @export
#'

compute_EC_fixed_number_windows <- function(table_daily_W = x,
                                            base_temperature = 10,
                                            method_GDD_calculation =
                                              c('method_b'),
                                            nb_windows_intervals = 8,
                                            ...) {
  
  # Calculation GDD
  table_daily_W$TMIN_GDD = table_daily_W$T2M_MIN
  table_daily_W$TMAX_GDD = table_daily_W$T2M_MAX
  
  if (method_GDD_calculation == 'method_b') {
    # Method b: when the minimum temperature T_min is below the T_base:
    # Any temperature below T_base is set to T_base before calculating the
    # average. https://en.wikipedia.org/wiki/Growing_degree-day
    
    table_daily_W$TMIN_GDD[table_daily_W$TMIN_GDD < base_temperature] <-
      base_temperature
    table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD < base_temperature] <-
      base_temperature
  }
  
  # The maximum temperature is usually capped at 30 °C for GDD calculation.
  
  table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD > 30] <- 30
  table_daily_W$TMEAN_GDD <-
    (table_daily_W$TMAX_GDD + table_daily_W$TMIN_GDD) / 2
  table_daily_W$GDD = table_daily_W$TMEAN_GDD - 10
  
  
  # Calculation day length
  
  table_daily_W$day_length <-
    daylength(lat = table_daily_W$LAT, day_of_year = table_daily_W$DOY)
  table_daily_W$PhotothermalTime <-
    table_daily_W$day_length * table_daily_W$GDD
  
  # Determine the width of each window based on the total number of windows
  # to use.
  duration_time_window_days <-
    floor(unique(table_daily_W$length.gs) / nb_windows_intervals)
  
  mean_TMIN <- zoo::rollapply(table_daily_W$T2M_MIN,
                              width = duration_time_window_days,
                              mean,
                              by = duration_time_window_days)
  
  mean_TMAX = zoo::rollapply(table_daily_W$T2M_MAX,
                             width = duration_time_window_days,
                             mean,
                             by = duration_time_window_days)
  
  
  mean_TMEAN = zoo::rollapply(table_daily_W$T2M,
                              width = duration_time_window_days,
                              mean,
                              by = duration_time_window_days)
  
  freq_TMAX_sup30 = zoo::rollapply(table_daily_W$T2M_MAX,
                                   width = duration_time_window_days,
                                   function(x) {
                                     length(which(x > 30)) / length(x)
                                   },
                                   by = duration_time_window_days)
  
  
  freq_TMAX_sup35 = zoo::rollapply(table_daily_W$T2M_MAX,
                                   width = duration_time_window_days,
                                   function(x) {
                                     length(which(x > 35)) / length(x)
                                   },
                                   by = duration_time_window_days)
  
  
  
  sum_GDD = zoo::rollapply(table_daily_W$GDD,
                           width = duration_time_window_days,
                           sum,
                           by = duration_time_window_days)
  
  sum_PTT = zoo::rollapply(table_daily_W$PhotothermalTime,
                           width = duration_time_window_days,
                           sum,
                           by = duration_time_window_days)
  
  sum_P = zoo::rollapply(table_daily_W$PRECTOT,
                         width = duration_time_window_days,
                         sum,
                         by = duration_time_window_days)
  
  freq_P_sup10 = zoo::rollapply(table_daily_W$PRECTOT,
                                width = duration_time_window_days,
                                function(x) {
                                  length(which(x > 10)) / length(x)
                                },
                                by = duration_time_window_days)
  
  sum_solar_radiation = zoo::rollapply(table_daily_W$ALLSKY_SFC_SW_DWN,
                                       width = duration_time_window_days,
                                       sum,
                                       by = duration_time_window_days)
  
  table_EC <-
    data.frame(
      mean_TMIN,
      mean_TMAX,
      mean_TMEAN,
      freq_TMAX_sup30,
      freq_TMAX_sup35,
      sum_GDD,
      sum_PTT,
      sum_P,
      freq_P_sup10,
      sum_solar_radiation
    )
  
  
  
  # Format for final EC table per environment
  # Each cell represents the value of the EC for this time window, e.g.
  # represents an EC on its own. Therefore, each cell should represent one
  # column.
  
  grid_tab <-
    as.data.frame(expand.grid(colnames(table_EC), row.names(table_EC)))
  grid_tab <- grid_tab[order(grid_tab$Var1), ]
  row.names(grid_tab) <- NULL
  
  table_EC_long <-
    data.frame(
      t(table_EC$mean_TMIN),
      t(table_EC$mean_TMAX),
      t(table_EC$mean_TMEAN),
      t(table_EC$freq_TMAX_sup30),
      t(table_EC$freq_TMAX_sup35),
      t(table_EC$sum_GDD),
      t(table_EC$sum_PTT),
      t(table_EC$sum_P),
      t(table_EC$freq_P_sup10),
      t(table_EC$sum_solar_radiation)
    )
  
  colnames(table_EC_long) <-
    paste0(grid_tab$Var1, '_', grid_tab$Var2)
  table_EC_long$IDenv <- unique(table_daily_W$IDenv)
  
  return(table_EC_long)
  
}








#' Compute the day length given the altitude and day of year.
#'
#' @param latitude \code{numeric} Latitude
#'
#' @param day_of_year \code{numeric} Day of year
#'
#' @return \code{numeric} Number of hours of daytime.
#'
#' @references  A model comparison for daylength as a function of latitude and
#'   day of year. Ecological Modelling, 80(1), 87-95. Forsythe, W. C., Rykiel Jr
#'   ,E. J., Stahl, R. S., Wu, H. I., & Schoolfield, R. M. (1995).
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export


daylength <- function(lat, day_of_year) {
  P = asin(.39795 * cos(.2163108 + 2 * atan(.9671396 * tan(
    .00860 * (day_of_year - 186)
  ))))
  
  D = 24 - (24 / pi) * acos((sin(0.8333 * pi / 180) + sin(lat * pi / 180) *
                               sin(P)) / (cos(lat * pi / 180) * cos(P)))
  return(D)
}