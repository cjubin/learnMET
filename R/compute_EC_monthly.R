#' Compute ECs on a monthly basis.
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
#'     \item sum_et0: number_total_fixed_windows columns, indicating the
#'     cumulative reference evapotranspiration over the respective day-window.
#'     \item freq_P_sup10: number_total_fixed_windows columns, indicating the
#'     frequency of days with total precipitation superior to 10 mm over the
#'     respective day-window.
#'     \item sum_solar_radiation: number_total_fixed_windows columns, indicating
#'     the accumulated incoming solar radiation over the respective day-window.
#'     \item mean_vapr_deficit: number_total_fixed_windows columns, indicating
#'     the mean vapour pressure deficit over the respective day-window.
#'     \item IDenv \code{character} ID of the environment (Location_Year)
#'    }
#'   @author Cathy C. Westhues \email{cathy.jubin@@hotmail.com}
#'   @export
#'


compute_EC_monthly <- function(table_daily_W = x,
                               base_temperature = 10,
                               method_GDD_calculation =
                                 c('method_b'),
                               capped_max_temperature = F,
                               max_temperature = 35,
                               ...) {
  checkmate::assert_names(
    colnames(table_daily_W),
    must.include  = c(
      'T2M_MIN',
      'T2M_MAX',
      'T2M',
      'daily_solar_radiation',
      'PRECTOTCORR'
    )
  )
  table_daily_W <-
    table_daily_W[order(as.Date(table_daily_W$YYYYMMDD)), ]
  
  if (!all(names(table_daily_W) %in% 'length.gs')) {
    table_daily_W$length.gs <- nrow(table_daily_W) - 1
  }
  
  
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
  
  # The maximum temperature can be capped at 30 °C for GDD calculation.
  
  if (capped_max_temperature) {
    table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD > max_temperature] <-
      max_temperature
    table_daily_W$TMIN_GDD[table_daily_W$TMIN_GDD > max_temperature] <-
      max_temperature
  }
  
  table_daily_W$TMEAN_GDD <-
    (table_daily_W$TMAX_GDD + table_daily_W$TMIN_GDD) / 2
  table_daily_W$GDD = table_daily_W$TMEAN_GDD - base_temperature
  
  
  # Calculation day length
  
  table_daily_W$day_length <-
    daylength(lat = table_daily_W$latitude, day_of_year = table_daily_W$DOY)
  table_daily_W$PhotothermalTime <-
    table_daily_W$day_length * table_daily_W$GDD
  
  
  
  mean_TMIN <- as.numeric(aggregate(T2M_MIN~MM+YEAR+IDenv, table_daily_W, mean)[,4])  
  
  mean_TMAX <- as.numeric(aggregate(T2M_MAX~MM+YEAR+IDenv, table_daily_W, mean)[,4])
  
  mean_TMEAN <- as.numeric(aggregate(T2M~MM+YEAR+IDenv, table_daily_W, mean)[,4])
  
  
  cumsum30_TMAX <- as.numeric(aggregate(T2M_MAX~MM+YEAR+IDenv, table_daily_W, function(x){
    sum(x[which(x > 30)])
  })[,4])
  
  
  freq_TMIN_inf_minus5 <- as.numeric(aggregate(T2M_MIN~MM+YEAR+IDenv, table_daily_W, function(x){
    length(which(x < (-5))) / length(x)
  })[,4])
  
  
  freq_TMAX_sup30 <- as.numeric(aggregate(T2M_MAX~MM+YEAR+IDenv, table_daily_W, function(x) {
    length(which(x > 30)) / length(x)
  })[,4])
  freq_TMAX_sup35 <- as.numeric(aggregate(T2M_MAX~MM+YEAR+IDenv, table_daily_W, function(x) {
    length(which(x > 35)) / length(x)
  })[,4])
  freq_TMAX_sup40 <- as.numeric(aggregate(T2M_MAX~MM+YEAR+IDenv, table_daily_W, function(x) {
    length(which(x > 40)) / length(x)
  })[,4])
  
  
  
  
  
  #sum_GDD = zoo::rollapply(table_daily_W$GDD,
  #                         width = duration_time_window_days,
  #                         sum,
  #                         by = duration_time_window_days)
  
  sum_PTT <- as.numeric(aggregate(table_daily_W$PhotothermalTime~MM+YEAR+IDenv, table_daily_W, sum)[,4])
  
  sum_P <- as.numeric(aggregate(table_daily_W$PRECTOTCORR~MM+YEAR+IDenv, table_daily_W, sum)[,4])
  
  
  
  if ("et0" %in% colnames(table_daily_W)) {
    sum_et0 <- as.numeric(aggregate(table_daily_W$et0~MM+YEAR+IDenv, table_daily_W, sum)[,4])      
  }
  
  mean_vapr_deficit <- as.numeric(aggregate(table_daily_W$vapr_deficit~MM+YEAR+IDenv, table_daily_W, mean)[,4])
  
  freq_P_sup10 <- as.numeric(aggregate(table_daily_W$PRECTOTCORR~MM+YEAR+IDenv, table_daily_W, function(x) {
    length(which(x > 10)) / length(x)
  })[,4])
  
  sum_solar_radiation <- as.numeric(aggregate(table_daily_W$daily_solar_radiation~MM+YEAR+IDenv, table_daily_W, sum)[,4])
  
  
  if ("et0" %in% colnames(table_daily_W)) {
    table_EC <-
      data.frame(
        mean_TMIN,
        mean_TMAX,
        mean_TMEAN,
        freq_TMAX_sup30,
        freq_TMAX_sup35,
        freq_TMAX_sup40,
        cumsum30_TMAX,
        #sum_GDD,
        sum_PTT,
        sum_P,
        sum_et0,
        freq_P_sup10,
        sum_solar_radiation,
        mean_vapr_deficit,
        freq_TMIN_inf_minus5
      )
  } else{
    table_EC <-
      data.frame(
        mean_TMIN,
        mean_TMAX,
        mean_TMEAN,
        freq_TMAX_sup30,
        freq_TMAX_sup35,
        freq_TMAX_sup40,
        cumsum30_TMAX,
        #sum_GDD,
        sum_PTT,
        sum_P,
        freq_P_sup10,
        sum_solar_radiation,
        mean_vapr_deficit,
        freq_TMIN_inf_minus5
      )
  }
  
  
  # Format for final EC table per environment
  # Each cell represents the value of the EC for this day-window, e.g.
  # represents an EC on its own. Therefore, each cell should represent one
  # column.
  
  grid_tab <-
    as.data.frame(expand.grid(colnames(table_EC), row.names(table_EC)))
  grid_tab <- grid_tab[order(grid_tab$Var1),]
  row.names(grid_tab) <- NULL
  
  
  if ("et0" %in% colnames(table_daily_W)) {
    table_EC_long <-
      data.frame(
        t(table_EC$mean_TMIN),
        t(table_EC$mean_TMAX),
        t(table_EC$mean_TMEAN),
        t(table_EC$freq_TMAX_sup30),
        t(table_EC$freq_TMAX_sup35),
        t(table_EC$freq_TMAX_sup40),
        t(table_EC$cumsum30_TMAX),
        #t(table_EC$sum_GDD),
        t(table_EC$sum_PTT),
        t(table_EC$sum_P),
        t(table_EC$sum_et0),
        t(table_EC$freq_P_sup10),
        t(table_EC$sum_solar_radiation),
        t(table_EC$mean_vapr_deficit),
        t(table_EC$freq_TMIN_inf_minus5)
      )
  }
  else{
    table_EC_long <-
      data.frame(
        t(table_EC$mean_TMIN),
        t(table_EC$mean_TMAX),
        t(table_EC$mean_TMEAN),
        t(table_EC$freq_TMAX_sup30),
        t(table_EC$freq_TMAX_sup35),
        t(table_EC$freq_TMAX_sup40),
        t(table_EC$cumsum30_TMAX),
        #t(table_EC$sum_GDD),
        t(table_EC$sum_PTT),
        t(table_EC$sum_P),
        t(table_EC$freq_P_sup10),
        t(table_EC$sum_solar_radiation),
        t(table_EC$mean_vapr_deficit),
        t(table_EC$freq_TMIN_inf_minus5)
      )
  }
  colnames(table_EC_long) <-
    paste0(grid_tab$Var1, '_', grid_tab$Var2)
  table_EC_long$IDenv <- unique(table_daily_W$IDenv)
  table_EC_long$year <- unique(table_daily_W$year)
  table_EC_long$location <- unique(table_daily_W$location)
  
  return(table_EC_long)
  
}


