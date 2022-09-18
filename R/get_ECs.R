#' Compute environmental covariates for each environment of the MET dataset.
#'
#' @description
#' This function enables to retrieve daily weather data from NASA POWER source
#' for each environment and derive environmental covariates over non-overlapping
#' time windows, which can be defined in various ways by the user.
#' The user can also provide own daily weather data, even for only part of the
#' total number of environments. For the remaining environments, weather data
#' will be retrieved using the NASA POWER query.
#'
#' @param info_environments a \code{data.frame} with the following columns:
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'     \item IDenv: \code{character} ID of the environment (location-year)
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD
#'  }

#'   * \strong{The data.frame should contain as many rows as Year x Location
#'   combinations which will be used in pheno_new.}
#'
#' @param raw_weather_data a \code{data.frame} which contains the following
#'   mandatory columns:
#'   \enumerate{
#'     \item longitude \code{numeric}
#'     \item latitude \code{numeric}
#'     \item year \code{numeric}
#'     \item location \code{character}
#'     \item YYYYMMDD \code{Date} Date of the daily observation
#'     \item IDenv \code{character} Environmt ID written Location_Year
#'     \item T2M \code{numeric} Average mean temperature (degree Celsius)
#'     \item T2M_MIN \code{numeric} Min. temperature (degree Celsius)
#'     \item T2M_MAX \code{numeric} Max. temperature (degree Celsius)
#'     \item PRECTOTCORR \code{numeric} Total daily precipitation (mm)
#'    }
#'   Additional weather data provided by user must be a subset of the following
#'   weather variable names (= next columns):
#'   (\strong{Any imputation step should be performed before providing
#'   this daily weather dataset to the package. }):
#'    \enumerate{
#'     \item RH2M \code{numeric} Daily mean relative humidity (%)
#'     \item RH2M_MIN \code{numeric} Daily minimum relative humidity (%)
#'     \item RH2M_MAX \code{numeric} Daily maximum relative humidity (%)
#'     \item daily_solar_radiation \code{numeric} daily solar radiation
#'     (MJ/m^2/day)
#'     \item T2MDEW \code{numeric} Dew Point (°C)
#'    }
#'    Default is `NULL`.
#'
#' @param method_ECs_intervals \code{character} A method among the four
#'   available:
#'   \enumerate{
#'     \item user_defined_intervals: if chosem the user must provide a data.frame
#'     intervals_growth_manual, described as parameter.
#'     \item GDD: day-intervals are determined based on growing degree days
#'     estimation, given a crop_model (argument must be provided).
#'     \item fixed_nb_windows_across_env: in each environment,
#'     the growing season is split into a number of windows equal to
#'     nb_windows_intervals.
#'     \item fixed_length_time_windows_across_env: in each environment,
#'     the growing season is divided into windows which always span the same
#'     length determined by the argument duration_time_window_days.
#'    }
#'   Default method is **fixed_nb_windows_across_env**.
#'
#'
#' @param fixed_length_time_windows_across_env \code{logical} indicates if the
#'   growing season lengths should be divided in non-overlapping time windows of
#'   fixed lengths (in days) across all environments.
#'   This implies that the total number of time windows, which need to be common
#'   across all environments, is determined by the shortest growing season
#'   included in the MET environments. This further implies that the total
#'   growing season may not be covered by the environmental predictors for
#'   the longest growing seasons.\cr
#'   Default is `TRUE`.
#
#'
#' @param fixed_nb_windows_across_env \code{logical} indicates if the
#'   growing season lengths should be divided in a fixed number of
#'   non-overlapping windows, which fully cover the growing season of each
#'   environment. This means that these time windows might not be of same length
#'   across environments (but always of same length in one environment). \cr
#'   Default is `FALSE`.
#'
#' @param nb_windows_intervals \code{numeric} Number of day-windows covering
#'   the growing season length (common number of day-windows across all
#'   environments). This argument is used if the default option for
#'   method_ECs_intervals is used ('fixed_nb_windows_across_env'). Default is 10.
#'
#' @param duration_time_window_days This argument is used only when the option
#'   'fixed_length_time_windows_across_env' is chosen. It determines the fixed
#'   number of days spanned within each window, across all environments.
#'   Default value is 10.
#'
#' @param base_temperature \code{numeric} It can be chosen by the user,
#'   to calculate GDD more accurately, based on the crop. Default value is 10
#'   degree Celsius. Base temperature will always be used by default.
#'
#' @param max_temperature \code{numeric} It can be chosen by the user,
#'   to calculate GDD by capping max temperature above this given threshold,
#'   based on the crop. Default value is 30 degree Celsius. By default, it is
#'   not used.
#'
#' @param crop_model \code{character} A crop_model among those implemented in
#'   [gdd_information()]. This argument is necessary only when the
#'   method_ECs_intervals called is "GDD". Default is NULL.
#'
#' @param intervals_growth_manual \code{data.frame} which is required only if
#'   the method_ECs_intervals chosen is "user_defined_intervals".
#'   * column 1: \code{numeric} year
#'   * column 2: \code{character} location
#'   * columns 3 and +: \code{numeric} Date (in Days after Planting) at which
#'   the crop enters a new growth stage in a given environment.
#'    "P" refers to the planting date and should contain 0 as value, "VE" to
#'    emergence, etc...
#'   \strong{Day 0 (Planting Date, denoted "P") should be in the third column.
#'   At least 4 columns should be in this data.frame. There is no need to
#'   indicate the column "Harvest" - already considered in the function.}
#'   An example of how this data.frame should be provided is given in
#'   [intervals_growth_manual_G2F].\cr
#'   Default is NULL.
#'
#' @param save_daily_weather_tables \code{logical} indicates whether the
#'   daily weather tables should be saved. Default is `TRUE`.
#'
#' @param path_data \code{character} Path of the folder where a
#'   RDS object will be created to save the daily weather tables if saved. (Do
#'   not use a Slash after the name of the last folder.)
#'
#' @param only_get_daily_data \code{logical} Only get daily weather data
#'
#' @param ... Arguments passed to the [compute_EC()] function.
#'
#' @return A \code{data.frame} object containing the weather-based environmental
#'   covariates.
#'
#' @references
#' \insertRef{sparks2018nasapower}{learnMET}
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#'
#' @export
#'
#'




get_ECs <-
  function(info_environments,
           raw_weather_data = NULL,
           method_ECs_intervals = "fixed_nb_windows_across_env",
           save_daily_weather_tables = T,
           path_data = NULL,
           crop_model = NULL,
           nb_windows_intervals = 10,
           duration_time_window_days = 10,
           base_temperature = 10,
           max_temperature = 30,
           capped_max_temperature = F,
           intervals_growth_manual = NULL,
           only_get_daily_data = FALSE,
           et0 = F,
           ...) {

    # Check the path_folder: create if does not exist



    if (!is.null(path_data)) {
      path_data <- file.path(path_data, "weather_data")
      if (!dir.exists(path_data)) {
        dir.create(path_data, recursive = T)
      }
    }



    if (is.null(info_environments$longitude) ||
      is.null(info_environments$latitude) ||
      is.na(info_environments$latitude) ||
      is.na(info_environments$longitude)) {
      stop("Longitude and latitude needed to impute ECs.\n")
    }

    if (is.null(info_environments$harvest.date) ||
      is.null(info_environments$planting.date) ||
      is.na(info_environments$harvest.date) ||
      is.na(info_environments$planting.date)) {
      stop("Planting and harvest dates needed to impute ECs (format Date, YYYY-MM-DD).\n")
    }

    # Checking that data are in the past to retrieve weather data

    assertive.datetimes::assert_all_are_in_past(x = info_environments$planting.date)
    assertive.datetimes::assert_all_are_in_past(x = info_environments$harvest.date)

    # Check if raw weather data for some environments are provided.
    # If yes, check which weather variables are provided.

    if (!is.null(raw_weather_data)) {
      if (is.null(path_data)) {
        stop(
          "Please indicate a path using path_to_save argument where the results from the QC can be stored."
        )
      }

      cat(
        paste(
          "Raw weather data are provided by the user and will be used",
          "to build environmental covariates.\n"
        )
      )

      raw_weather_data$IDenv <-
        paste0(raw_weather_data$location, "_", raw_weather_data$year)
      raw_weather_data$DOY <- as.integer(lubridate::yday(raw_weather_data$YYYYMMDD))
      raw_weather_data <-
        qc_raw_weather_data(
          daily_weather_data = raw_weather_data,
          info_environments = info_environments,
          et0 = et0,
          path_flagged_values = path_data
        )

      variables_raw_data <-
        colnames(raw_weather_data)
      list_envs_to_retrieve_all_data <-
        info_environments$IDenv[which(info_environments$IDenv %notin% raw_weather_data$IDenv)]
    } else {
      variables_raw_data <- NULL
      list_envs_to_retrieve_all_data <-
        unique(info_environments$IDenv)
    }



    ############################################################################
    # Obtain daily "AG" community daily weather information for each environment
    # using nasapower R package
    ############################################################################
    if (checkmate::test_character(list_envs_to_retrieve_all_data, all.missing = F)) {
      # Check that the data have not been downloaded before (via learnMET) and saved as RDS file
      # Also check that in that case, the planting and harvest dates are matching those presently used
      if (file.exists(file.path(
        path_data,
        "daily_weather_tables_nasapower.RDS"
      ))) {
        requested_data <- as.data.frame(readRDS(
          file.path(
            path_data,
            "daily_weather_tables_nasapower.RDS"
          )
        ))

        no_missing_env_previous <-
          checkmate::testNames(unique(requested_data$IDenv),
            must.include = unique(info_environments$IDenv)
          )
        if (!no_missing_env_previous) {
          list_envs_to_retrieve_completely <-
            unique(info_environments$IDenv)[which(unique(info_environments$IDenv) %notin% unique(requested_data$IDenv))]
        } else {
          list_envs_to_retrieve_completely <- NULL
        }

        # In case the previous run has more environments than the present analysis:
        requested_data <-
          requested_data[requested_data$IDenv %in% info_environments$IDenv, ]

        are_previous_data_OK <- vector()
        n <- 1
        for (j in unique(requested_data$IDenv)) {
          int <-
            lubridate::interval(info_environments[info_environments$IDenv == j, "planting.date"], info_environments[info_environments$IDenv ==
              j, "harvest.date"])

          are_previous_data_OK[n] <-
            ifelse(!all(lubridate::`%within%`(requested_data[requested_data$IDenv == j, "YYYYMMDD"],
              int)), FALSE, TRUE)
          n <- n + 1
        }

        env_to_keep <-
          unique(requested_data$IDenv)[which(are_previous_data_OK == TRUE)]

        list_envs_to_redownload <-
          unique(requested_data$IDenv)[which(are_previous_data_OK == FALSE)]

        data_previous_run <-
          split(requested_data[requested_data$IDenv %in% env_to_keep, ], f = requested_data$IDenv)
        list_envs_to_retrieve_all_data <-
          c(
            list_envs_to_retrieve_completely,
            list_envs_to_redownload
          )
        cat(
          "Daily weather tables have been downloaded from NASA POWER for the required environments in a previous run, and are matching the environments ID/planting and harvest dates used in this analysis.\n These data will be used. \n"
        )
      } else {
        list_envs_to_retrieve_completely <- NULL
        list_envs_to_redownload <- NULL
        data_previous_run <- NULL
      }

      # If required, retrieving and saving soil data for all environments

      # requested_data_soil <- tryCatch({
      #  get_soil_per_env(
      #    environment = environment,
      #    info_environments = info_environments
      #  )
      #
      # },
      # error = function(e)
      #  return(NULL),
      # warning = function(w)
      #  return(NULL))

      # If we do not have any previously saved data, or if some environments
      # are missing/incorrect (= dates in the growing season missing) for the
      # new analysis:
      if (!file.exists(file.path(
        path_data,
        "daily_weather_tables_nasapower.RDS"
      )) |
        !is.null(list_envs_to_redownload) |
        !is.null(list_envs_to_retrieve_completely)) {
        has_unsuccessful_requests <- TRUE
        counter <- 1

        list_envs_loop <- list_envs_to_retrieve_all_data
        # This is an empty list to which all requested data will be assigned.
        requested_data <-
          vector(
            mode = "list",
            length = length(list_envs_loop)
          )
        names(requested_data) <- list_envs_loop

        # Issues with the NASAPOWER query: it sometimes fail --> use of tryCath
        # and while procedure to ensure weather data for each envrionment
        # are retrieved.
        while (has_unsuccessful_requests) {
          res_w_daily_all <-
            lapply(
              list_envs_loop,
              function(environment, ...) {
                requested_data <- tryCatch(
                  {
                    get_daily_tables_per_env(
                      environment = environment,
                      info_environments = info_environments,
                      path_data = path_data,
                      et0 = et0,
                      ...
                    )
                  },
                  error = function(e) {
                    return(NULL)
                  },
                  warning = function(w) {
                    return(NULL)
                  }
                )

                return(requested_data)
              }
            )
          names(res_w_daily_all) <- list_envs_loop
          unsuccessful_request_bool <- vapply(res_w_daily_all,
            FUN = is.null,
            FUN.VALUE = logical(1)
          )

          failed_requests <-
            list_envs_loop[unsuccessful_request_bool]
          good_requests <-
            list_envs_loop[!unsuccessful_request_bool]

          list_envs_loop <- failed_requests


          requested_data[good_requests] <-
            res_w_daily_all[good_requests]

          counter <- counter + 1

          if (counter == 15) {
            stop("At least one request failed fifteen times.", call. = FALSE)
          }

          has_unsuccessful_requests <- any(unsuccessful_request_bool)
        }

        if (!is.null(data_previous_run)) {
          requested_data <- c(data_previous_run, requested_data)
        }

        cat(
          "Daily weather tables downloaded from NASA POWER for the required",
          "environments!\n"
        )
      }

      # Save daily weather data used to compute ECs

      if (save_daily_weather_tables & !is.null(path_data)) {
        saveRDS(
          data.table::rbindlist(requested_data),
          file.path(
            path_data,
            "daily_weather_tables_nasapower.RDS"
          )
        )
      }
      climate_data_retrieved <- TRUE
    } else {
      climate_data_retrieved <- FALSE
    }

    #######################################################################
    ## Pre-processing to merge user data + NASA data
    #######################################################################


    if (is.null(raw_weather_data)) {
      # All weather date were retrieved from NASAPOWER data source
      weather_data_list <- requested_data
    }

    if (!is.null(raw_weather_data)) {
      if (length(list_envs_to_retrieve_all_data) > 1) {
        # Missing weather information for some environments is binded to
        # weather information provided by the user for some environments.

        raw_weather_data$IDenv <- as.factor(raw_weather_data$IDenv)
        weather_data_list <-
          split(raw_weather_data, raw_weather_data$IDenv)
        weather_data_list <-
          append(weather_data_list, requested_data)
      }

      if (length(list_envs_to_retrieve_all_data) == 0) {
        # All weather date were provided by the user as daily weather data tables.

        raw_weather_data$IDenv <- as.factor(raw_weather_data$IDenv)
        weather_data_list <-
          split(raw_weather_data, raw_weather_data$IDenv)
      }
    }

    if (only_get_daily_data) {
      return(weather_data_list)
    }

    #############################################
    # Derivation of EC based on selected method #
    #############################################


    if (method_ECs_intervals == "user_defined_intervals") {
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_user_defined_intervals(
              table_daily_W = x,
              intervals_growth_manual = intervals_growth_manual,
              base_temperature = base_temperature,
              max_temperature = max_temperature,
              capped_max_temperature = capped_max_temperature,
              ...
            )
          }
        )


      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs <-
        merged_ECs[, c("IDenv", "year", "location", colnames(merged_ECs)[colnames(merged_ECs) %notin%
          c("IDenv", "year", "location")])]
    }

    if (method_ECs_intervals == "GDD") {
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_gdd(
              table_daily_W = x,
              crop_model = crop_model,
              capped_max_temperature = capped_max_temperature,
              ...
            )
          }
        )


      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs <-
        merged_ECs[, c("IDenv", "year", "location", colnames(merged_ECs)[colnames(merged_ECs) %notin%
          c("IDenv", "year", "location")])]
    }

    if (method_ECs_intervals == "fixed_length_time_windows_across_env") {
      # Each EC is computed over a fixed certain number of days, given by the
      # parameter "duration_time_window_days".
      # The maximum number of time windows (e.g. the total number of ECs)
      # is determined by the shortest growing season across all environments.

      length_minimum_gs <- min(vapply(weather_data_list, function(x) {unique(as.numeric(x[, "length.gs"]))}, numeric(1)))

      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_fixed_length_window(
              table_daily_W = x,
              length_minimum_gs = length_minimum_gs,
              duration_time_window_days = duration_time_window_days,
              base_temperature = base_temperature,
              max_temperature = max_temperature,
              capped_max_temperature = capped_max_temperature,
              ...
            )
          }
        )



      merged_ECs <- do.call("rbind", ECs_all_envs)
      merged_ECs <-
        merged_ECs[, c("IDenv", "year", "location", colnames(merged_ECs)[colnames(merged_ECs) %notin%
          c("IDenv", "year", "location")])]
    }


    if (method_ECs_intervals == "fixed_nb_windows_across_env") {
      ECs_all_envs <-
        lapply(
          weather_data_list,
          FUN = function(x, ...) {
            compute_EC_fixed_number_windows(
              table_daily_W = x,
              nb_windows_intervals = nb_windows_intervals,
              base_temperature = base_temperature,
              max_temperature = max_temperature,
              capped_max_temperature = capped_max_temperature,
              ...
            )
          }
        )


      merged_ECs <- do.call("rbind", ECs_all_envs)

      merged_ECs <-
        merged_ECs[, c("IDenv", "year", "location", colnames(merged_ECs)[colnames(merged_ECs) %notin%
          c("IDenv", "year", "location")])]
    }

    merged_ECs <- list(
      "ECs" = merged_ECs,
      "climate_data_retrieved" = climate_data_retrieved
    )



    return(merged_ECs)
  }
