compute_EC_wofost <- function(table_daily_W,
                           crop_model,
                           method_GDD_calculation =
                             c('method_b'),
                           ...) {
  
  
  
  ## Define days for which a new stage is reached in terms of GDD
  
  new_stage_reached <-unlist(lapply(table_gdd$GDD,FUN = function(x){min(which(table_daily_W$cumGDD>x))}))
  new_stage_reached <- c(0,new_stage_reached)
  
  table_daily_W$interval=cut(seq_len(nrow(table_daily_W)), breaks = new_stage_reached, include.lowest = TRUE, right = FALSE)
  
  
  mean_TMIN <- unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) mean(x$T2M_MIN)))
  
  mean_TMAX = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) mean(x$T2M_MAX)))
  
  mean_TMEAN = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) mean(x$T2M)))
  
  freq_TMAX_sup30 = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) {
    length(which(x$T2M_MAX > 30)) / length(x$T2M_MAX)
  }))
  
  
  freq_TMAX_sup35 = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) {
    length(which(x$T2M_MAX > 35)) / length(x$T2M_MAX)
  }))
  
  
  sum_PTT = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) sum(x$PhotothermalTime)))
  
  
  sum_P =  unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) sum(x$PRECTOT)))
  
  freq_P_sup10 = unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) {
    length(which(x$PRECTOT > 30)) / length(x$PRECTOT)
  }))
  
  sum_solar_radiation =  unlist(lapply(split(table_daily_W, f = table_daily_W$interval), FUN = function(x) sum(x$ALLSKY_SFC_SW_DWN)))
  
  
  table_EC <-
    data.frame(
      mean_TMIN,
      mean_TMAX,
      mean_TMEAN,
      freq_TMAX_sup30,
      freq_TMAX_sup35,
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