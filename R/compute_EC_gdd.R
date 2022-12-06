#' Compute ECs based on growth stages which are estimated based on accumulated
#' GDD in each environment.
#'
#' @description
#' This function enables to retrieve daily weather data for each
#' environment and derive environmental covariates over non-overlapping time
#' windows, which can be defined in various ways by the user.
#'
#' @param table_daily_W \code{data.frame} Object returned by the function
#'   [get_daily_tables_per_env()]
#'
#' @param crop_model \code{character} Name of the crop model used to estimate
#'   the times of the crop stages based on temperature sum accumulation.
#'   Current options are `maizehybrid1700` and `hardwheatUS`. Growing degree
#'   days are utilized to delineate maize phenology.
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
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export

compute_EC_gdd <- function(table_daily_W,
                           crop_model = NULL,
                           method_GDD_calculation =
                             c("method_b"),
                           capped_max_temperature = F,
                           ...) {
  checkmate::assert_character(crop_model)
  checkmate::assert_names(
    colnames(table_daily_W),
    must.include  = c("T2M_MIN",
                      "T2M_MAX",
                      "T2M",
                      "daily_solar_radiation")
  )
  
  table_daily_W <-
    table_daily_W[order(as.Date(table_daily_W$YYYYMMDD)), ]
  
  
  table_gdd <- gdd_information(crop_model = crop_model)[[1]]
  base_temperature <- gdd_information(crop_model = crop_model)[[2]]
  
  if (crop_model %in% c("barley_hawn")) {
    max_temperature1 <- gdd_information(crop_model = crop_model)[[3]]
    max_temperature2 <-
      gdd_information(crop_model = crop_model)[[4]]
    stage_change_max_temp <-
      gdd_information(crop_model = crop_model)[[5]]
    
  }
  else{
    max_temperature <- gdd_information(crop_model = crop_model)[[3]]
  }
  
  
  
  # Calculation GDD
  table_daily_W$TMIN_GDD <- table_daily_W$T2M_MIN
  table_daily_W$TMAX_GDD <- table_daily_W$T2M_MAX
  
  
  if (method_GDD_calculation == "method_b") {
    # Method b: when the minimum temperature T_min is below the T_base:
    # Any temperature below T_base is set to T_base before calculating the
    # average. https://en.wikipedia.org/wiki/Growing_degree-day
    
    table_daily_W$TMIN_GDD[table_daily_W$TMIN_GDD < base_temperature] <-
      base_temperature
    table_daily_W$TMAX_GDD[table_daily_W$TMAX_GDD < base_temperature] <-
      base_temperature
    
  }
  
  if (crop_model %in% c("barley_hawn")) {
    table_daily_W$TMAX_GDD_before_2 <- table_daily_W$TMAX_GDD
    table_daily_W$TMAX_GDD_after_2 <- table_daily_W$TMAX_GDD
  }
  
  # The maximum temperature can be capped for GDD calculation.
  
  if (capped_max_temperature) {
    if (crop_model %in% c("barley_hawn")) {
      table_daily_W$TMAX_GDD_before_2[which(table_daily_W$TMAX_GDD_before_2 > max_temperature1)] <-
        max_temperature1
      table_daily_W$TMAX_GDD_after_2[which(table_daily_W$TMAX_GDD_after_2 > max_temperature2)] <-
        max_temperature2
      
    } else{
      table_daily_W$TMAX_GDD[which(table_daily_W$TMAX_GDD > max_temperature)] <-
        max_temperature
      table_daily_W$TMIN_GDD[which(table_daily_W$TMIN_GDD > max_temperature)] <-
        max_temperature
    }
    
    
  }
  
  if (crop_model %notin% c("barley_hawn")) {
    table_daily_W$TMEAN_GDD <-
      (table_daily_W$TMAX_GDD + table_daily_W$TMIN_GDD) / 2
    table_daily_W$GDD = table_daily_W$TMEAN_GDD - base_temperature
    
  } else{
    table_daily_W$TMEAN_GDD_before_2 <-
      (table_daily_W$TMAX_GDD_before_2 + table_daily_W$TMIN_GDD) / 2
    table_daily_W$GDD_before_2 = table_daily_W$TMEAN_GDD_before_2 - base_temperature
    
    table_daily_W$TMEAN_GDD_after_2 <-
      (table_daily_W$TMAX_GDD_after_2 + table_daily_W$TMIN_GDD) / 2
    table_daily_W$GDD_after_2 = table_daily_W$TMEAN_GDD_after_2 - base_temperature
    
    change_stage <-
      as.numeric(table_gdd[which(table_gdd$Stage == stage_change_max_temp), "GDD"])
    
    table_daily_W$cumGDD_before_2 <-
      cumsum(table_daily_W$GDD_before_2)
    
    table_daily_W$GDD <- NA
    table_daily_W$GDD[which(table_daily_W$cumGDD_before_2 < change_stage)] <-
      table_daily_W$GDD_before_2[which(table_daily_W$cumGDD_before_2 < change_stage)]
    table_daily_W$GDD[which(table_daily_W$cumGDD_before_2 > change_stage)] <-
      table_daily_W$GDD_after_2[which(table_daily_W$cumGDD_before_2 > change_stage)]
    
  }
  
  
  if (method_GDD_calculation == "method_a") {
    table_daily_W$GDD[table_daily_W$GDD < 0] <- 0
  }
  
  
  # Calculation day length
  if ("SG_DAY_HOUR_AVG" %in% names(table_daily_W)) {
    table_daily_W$day_length <- table_daily_W$SG_DAY_HOUR_AVG
  }
  else{
    table_daily_W$day_length <-
      daylength(lat = table_daily_W$latitude, day_of_year = table_daily_W$DOY)
  }
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
    print(paste0(
      "GDDs missing for the environment",
      unique(table_daily_W$IDenv)
    ))
    new_stage_reached <-
      new_stage_reached[-which(new_stage_reached == Inf)]
    new_stage_reached <- c(new_stage_reached, nrow(table_daily_W))
  }
  
  new_stage_reached <- c(0, new_stage_reached, nrow(table_daily_W))
  
  table_daily_W$interval = cut(
    seq_len(nrow(table_daily_W)),
    breaks = new_stage_reached,
    include.lowest = TRUE,
    right = FALSE
  )
  
  
  intervals_growth <- c(0, table_gdd$Stage, "Harvest")
  levels(table_daily_W$interval) <-
    paste(intervals_growth[1:(length(intervals_growth) - 1)], intervals_growth[2:(length(intervals_growth))], sep = "-")
  
  
  if ("T2M_MIN" %in% names(table_daily_W)) {
    mean_TMIN <-
      unlist(lapply(
        split(table_daily_W, f = table_daily_W$interval),
        FUN = function(x)
          mean(x$T2M_MIN)
      ))
    freq_TMIN_inf_minus5 = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x) {
        length(which(x$T2M_MIN < (-5))) / length(x$T2M_MIN)
      }
    ))
    
    
    
  }
  
  if ("T2M_MAX" %in% names(table_daily_W)) {
    mean_TMAX = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        mean(x$T2M_MAX)
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
    freq_TMAX_sup40 = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x) {
        length(which(x$T2M_MAX > 40)) / length(x$T2M_MAX)
      }
    ))
    
    cumsum30_TMAX = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x) {
        sum(x[which(x$T2M_MAX > 30), "T2M_MAX"])
      }
    ))
  }
  if ("T2M" %in% names(table_daily_W)) {
    mean_TMEAN = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        mean(x$T2M)
    ))
  }
  
  if ("PhotothermalTime" %in% names(table_daily_W)) {
    sum_PTT = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        sum(x$PhotothermalTime, na.rm = T)
    ))
  }
  if ("et0" %in% colnames(table_daily_W)) {
    sum_et0 = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        sum(x$et0, na.rm = T)
    ))
  }
  
  if ("PRECTOTCORR" %in% names(table_daily_W)) {
    sum_P =  unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        sum(x$PRECTOTCORR, na.rm = T)
    ))
    
    freq_P_sup10 = unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x) {
        length(which(x$PRECTOTCORR > 10)) / length(x$PRECTOTCORR)
      }
    ))
  }
  
  if ("daily_solar_radiation" %in% names(table_daily_W)) {
    sum_solar_radiation =  unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        sum(x$daily_solar_radiation, na.rm = T)
    ))
  }
  
  if ("vapr_deficit" %in% names(table_daily_W)) {
    mean_vapr_deficit =  unlist(lapply(
      split(table_daily_W, f = table_daily_W$interval),
      FUN = function(x)
        mean(x$vapr_deficit, na.rm = T)
    ))
  }
  toMatch <- c("mean|freq_|sum_|cumsum")
  matches <- ls(pattern = toMatch)
  table_EC <-
    as.data.frame(`row.names<-`(do.call(cbind, mget(matches)), NULL))
  View(table_EC)
  
  row.names(table_EC) <- 1:nrow(table_EC)
  
  
  # Format for final EC table per environment
  # Each cell represents the value of the EC for this day-window, e.g.
  # represents an EC on its own. Therefore, each cell should represent one
  # column.
  
  table_EC_long <- unlist(table_EC)
  
  table_EC_long$IDenv <- unique(table_daily_W$IDenv)
  table_EC_long$year <- unique(table_daily_W$year)
  table_EC_long$location <- unique(table_daily_W$location)
  
  
  return(table_EC_long)
  
}