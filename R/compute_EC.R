compute_EC <- function(table_daily_W,
                       duration_time_window_days,
                       base_temperature,
                       method_GDD_calculation,
                       number_total_fixed_windows) {
  
  # Calculation GDD
  table_daily_W$TMIN_GDD = table_daily_W$T2M_MIN
  table_daily_W$TMAX_GDD = table_daily_W$T2M_MAX
  
  if (method_GDD_calculation == 'method_b') {
    # Method b: when the minimum temperature T_min is below the T_base:
    # Any temperature below T_base is set to T_base before calculating the average
    # https://en.wikipedia.org/wiki/Growing_degree-day
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
  # Each cell represents the value of the EC for this time window, e.g. represents
  # an EC on its own. Therefore, each cell should represent one column.
  
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
  
  colnames(table_EC_long) <- paste0(grid_tab$Var1, '_', grid_tab$Var2)
  table_EC_long$IDenv <- unique(table_daily_W$IDenv)
  
  return(table_EC_long)
  
}


# Forsythe, W. C., Rykiel Jr, E. J., Stahl, R. S., Wu, H. I., & Schoolfield, R. M. (1995). 
# A model comparison for daylength as a function of latitude and day of year. Ecological Modelling, 80(1), 87-95.

daylength = function(lat, day_of_year) {
  P = asin(.39795 * cos(.2163108 + 2 * atan(.9671396 * tan(
    .00860 * (day_of_year - 186)
  ))))
  
  D = 24 - (24 / pi) * acos((sin(0.8333 * pi / 180) + sin(lat * pi / 180) *
                               sin(P)) / (cos(lat * pi / 180) * cos(P)))
  return(D)
}