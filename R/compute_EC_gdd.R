#' Compute ECs based on growth requirements estimated from accumulated GDD.
#'
#' @description
#' This function enables to retrieve daily weather data for each
#' environment and derive environmental covariates over non-overlapping time
#' windows, which can be defined in various ways by the user.
#' @param table_daily_W
#' @param crop_model
#' @param method_GDD_calculation
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export

compute_EC_gdd <- function(table_daily_W,
                           crop_model,
                           method_GDD_calculation =
                             c('method_b'),
                           ...) {
  table_gdd <- gdd_information(crop_model = crop_model)[[1]]
  base_temperature <- gdd_information(crop_model = crop_model)[[2]]
  
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
  
  table_daily_W$cumGDD <- cumsum(table_daily_W$GDD)
  
  ## Define days for which a new stage is reached in terms of GDD
  
  new_stage_reached <-
    unlist(lapply(
      table_gdd$GDD,
      FUN = function(x) {
        min(which(table_daily_W$cumGDD > x))
      }
    ))
  
  if (Inf %in% new_stage_reached) {
    new_stage_reached <-
      new_stage_reached[-which(new_stage_reached == Inf)]
    new_stage_reached <- c(new_stage_reached,nrow(table_daily_W))
  }
  
  new_stage_reached <- c(0, new_stage_reached)
  
  table_daily_W$interval = cut(
    seq_len(nrow(table_daily_W)),
    breaks = new_stage_reached,
    include.lowest = TRUE,
    right = FALSE
  )
  
  
  mean_TMIN <-
    unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        mean(x$T2M_MIN)
    ))
  
  mean_TMAX = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x)
      mean(x$T2M_MAX)
  ))
  
  mean_TMEAN = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x)
      mean(x$T2M)
  ))
  
  freq_TMAX_sup30 = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x) {
      length(which(x$T2M_MAX > 30)) / length(x$T2M_MAX)
    }
  ))
  
  
  freq_TMAX_sup35 = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x) {
      length(which(x$T2M_MAX > 35)) / length(x$T2M_MAX)
    }
  ))
  
  
  sum_PTT = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x)
      sum(x$PhotothermalTime)
  ))
  
  
  sum_P =  unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x)
      sum(x$PRECTOT)
  ))
  
  freq_P_sup10 = unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x) {
      length(which(x$PRECTOT > 30)) / length(x$PRECTOT)
    }
  ))
  
  sum_solar_radiation =  unlist(lapply(
    split(table_daily_W, f = table_daily_W$interval),
    FUN = function(x)
      sum(x$ALLSKY_SFC_SW_DWN)
  ))
  
  
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
    table_EC <- table_EC[1:number_total_fixed_windows,]
  }
  
  # Format for final EC table per environment
  # Each cell represents the value of the EC for this time window, e.g.
  # represents an EC on its own. Therefore, each cell should represent one
  # column.
  
  grid_tab <-
    as.data.frame(expand.grid(colnames(table_EC), row.names(table_EC)))
  grid_tab <- grid_tab[order(grid_tab$Var1),]
  row.names(grid_tab) <- NULL
  
  table_EC_long <-
    data.frame(
      t(table_EC$mean_TMIN),
      t(table_EC$mean_TMAX),
      t(table_EC$mean_TMEAN),
      t(table_EC$freq_TMAX_sup30),
      t(table_EC$freq_TMAX_sup35),
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