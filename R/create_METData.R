#' Create a multi-environment trials data object
#'
#' @description
#' This function combines all types of data sources (genotypic, phenotypic,
#' information about the environments, environmental data if available...)
#' in a single data object of class \code{METData}.
#'
#' @name create_METData
#'
#' @param geno \code{numeric} genotype values stored in a \code{matrix} or
#'   \code{data.frame} which contains the geno_ID as row.names and markers as
#'   columns.
#'
#' @param map \code{data.frame} object with 3 columns.
#'   \enumerate{
#'   \item marker \code{character} with marker names
#'   \item chr \code{numeric} with chromosome number
#'   \item pos \code{numeric} with marker position.
#'   }
#'   \strong{Map object not mandatory}.
#'
#' @param pheno \code{data.frame} object with at least 4 columns.
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'   }
#'   From the fourth column on: each column is \code{numeric} and contains
#'   phenotypic values for a phenotypic trait observed in a combination
#'   Year x Location. Names of the traits should be provided as column names.
#'   \cr
#'   * \strong{The geno_ID must be a subset of the row.names in the geno object.
#'   }
#'
#' @param info_environments \code{data.frame} object with at least
#'   the 4 following columns. \cr
#'    \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'  }
#'  The two next columns are required only if weather data should be
#'  retrieved from NASA POWER data using the argument `compute_climatic_EC` set
#'  to TRUE, or if raw weather data are provided:
#'  \enumerate{
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD 
#'     \item elevation: (optional) \code{numeric} \cr
#'   }
#'   * \strong{The data.frame should contain as many rows as Year x Location
#'   combinations which will be used in pheno_new.}
#'
#' @param climate_variables \code{data.frame} can be let as NULL by user, if no
#'   climate variables provided as input. Otherwise, a \code{data.frame} should
#'   be provided.
#'   \strong{The data.frame should contain as many rows as the `info_environments`
#'    \code{data.frame}.} \cr
#'   Columns should be:
#'   \enumerate{
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'   }
#'   Columns 3 and + should be numeric and contain the climate (weather-based)
#'   covariates.\cr
#'
#'   * \strong{If climate_variables is provided,`compute_climatic_ECs`should be
#'   set to `FALSE`.}
#'
#' @param soil_variables \code{data.frame} can be let as NULL by user, if no
#'   soil variables provided as input. Otherwise, a \code{data.frame} should
#'   be provided.
#'   \strong{The data.frame should contain as many rows as the `info_environments`
#'    \code{data.frame}.} \cr
#'   Columns should be:
#'   \enumerate{
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'   }
#'   Columns 3 and + should be numeric and contain the soil-based environmental
#'   covariates.\cr
#'
#' @param raw_weather_data \code{data.frame} can be let as NULL by user, if no
#'   daily weather datasets are available. If else, required columns should be
#'   provided like this (colnames should be respected):
#'   \enumerate{
#'     \item longitude \code{numeric}
#'     \item latitude \code{numeric}
#'     \item year \code{numeric}
#'     \item location \code{character}
#'     \item YYYYMMDD \code{Date}
#'   }
#'   Available weather data provided by user must be a subset of the following
#'   weather variable names. Colnames must be given as following:
#'   \enumerate{
#'     \item T2M \code{numeric} Daily mean temperature (°C)
#'     \item T2M_MIN \code{numeric} Daily minimum temperature (°C)
#'     \item T2M_MAX \code{numeric} Daily maximum temperature (°C)
#'     \item PRECTOTCORR \code{numeric} Daily total precipitation (mm)
#'     \item RH2M \code{numeric} Daily mean relative humidity (%)
#'     \item RH2M_MIN \code{numeric} Daily minimum relative humidity (%)
#'     \item RH2M_MAX \code{numeric} Daily maximum relative humidity (%)
#'     \item daily_solar_radiation \code{numeric} daily solar radiation
#'     (MJ/m^2/day)
#'     \item top_atmosphere_insolation \code{numeric} Top-of-atmosphere
#'     Insolation (MJ/m^2/day)
#'     \item T2MDEW \code{numeric} Dew Point (°C)
#'    }
#'
#'   \strong{It is not required that weather data for ALL environments are
#'   provided by the user. If weather data for some environments are missing,
#'   they will be retrieved by the NASA }
#'
#'
#' @param compute_climatic_ECs \code{logical} indicates if climatic covariates
#'   should be computed with the function. Default
#'   is `FALSE`. \cr
#'   \strong{Set compute_climatic_ECs = `TRUE` if user wants to use weather data
#'   from NASA POWER data OR if raw weather data are available and should be
#'   used (also possible to provide field weather data for only some
#'   environments; weather data for other environments present in the dataset 
#'   will be retrieved using the NASA POWER query.}
#'
#' @param path_to_save Path where daily weather data (if retrieved) and plots 
#'   based on k-means clustering are saved.
#'   
#' @param as_test_set If using a prediction set (i.e. no phenotypic values
#'   for the new data to predict), should be set to TRUE. Default is FALSE.
#'
#' @return A formatted \code{list} of class \code{METData} which contains the
#'   following elements:
#'
#' * **geno**: \code{matrix} with genotype values of phenotyped individuals.
#'
#' * **map**: \code{data.frame} with genetic map.
#'
#' * **pheno**: \code{data.frame} with phenotypic trait values.
#'
#' * **compute_EC_by_geno**: \code{logical} indicates if environmental
#'   covariates were required to be retrieved via the package by the user.
#'
#' * **env_data**: \code{data.frame} with the environmental covariates per
#'   environment
#'
#' * **list_climatic_predictors**: \code{character} with the names of the climatic predictor variables
#'
#' * **list_soil_predictors**: \code{character} with the names of the soil-based predictor variables
#'
#' * **info_environments**: \code{data.frame} contains basic information on
#'   each environment.
#'
#' * **ECs_computed**: \code{logical} subelement added in the output
#'   to indicate if the function [get_ECs()] was run within the pipeline.
#'
#' * **climate_data_retrieved**: \code{logical} subelement added in the output
#'   to indicate if NASAPOWER data were retrieved within the pipeline.
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' @examples
#'
#' data(geno_G2F)
#' data(pheno_G2F)
#' data(map_G2F)
#' data(info_environments_G2F)
#' data(soil_G2F)
#' # Create METData and get climate variables from NASAPOWER data & use soil variables
#' METdata_G2F <- create_METData(geno=geno_G2F,pheno=pheno_G2F,map=map_G2F,climate_variables = NULL,compute_climatic_ECs = TRUE,info_environments = info_environments_G2F,soil_variables=soil_G2F, path_to_save = "~/g2f_data")
#'
#' data(geno_indica)
#' data(map_indica)
#' data(pheno_indica)
#' data(info_environments_indica)
#' data(climate_variables_indica)
#' METdata_indica <- create_METData(geno=geno_indica,pheno=pheno_indica,climate_variables = climate_variables_indica,compute_climatic_ECs = FALSE,info_environments = info_environments_indica,map = map_indica, path_to_save = "~/indica")
#'
#' data(geno_japonica)
#' data(map_japonica)
#' data(pheno_japonica)
#' data(info_environments_japonica)
#' data(climate_variables_japonica)
#' METdata_japonica <- create_METData(geno=geno_japonica,pheno=pheno_japonica,climate_variables = climate_variables_japonica,compute_climatic_ECs = FALSE,info_environments = info_environments_japonica,map = map_japonica, path_to_save = "~/japonica")

new_create_METData <-
  function(geno = NULL,
           map = NULL,
           pheno = NULL,
           info_environments = NULL,
           raw_weather_data = NULL,
           climate_variables = NULL,
           soil_variables = NULL,
           compute_climatic_ECs = FALSE,
           path_to_save = NULL,
           as_test_set = FALSE,
           get_public_soil_data = FALSE,
           ...) {
    # check if one object is missing / appropriate classes
    
    # If geno provided as data.frame --> convert to matrix after check on residual missing values
    
    if (is.data.frame(geno)) {
      checkmate::assert_data_frame(geno, any.missing = F, types = "numeric")
      geno <- as.matrix(geno)
    }
    
    
    checkmate::assert_matrix(geno, any.missing = F, mode = "numeric")
    
    if (!as_test_set) {
      checkmate::assert_data_frame(pheno, all.missing = F, min.cols = 4)
      
      checkmate::assert_names(names(pheno),
                              must.include = c("geno_ID",
                                               "year",
                                               "location"))
    } else{
      checkmate::assert_data_frame(pheno, all.missing = F, min.cols = 3)
    }
    
    checkmate::assert_data_frame(map, null.ok = T)
    
    if (!is.null(map)) {
      checkmate::assert_names(names(map),
                              must.include = c("marker",
                                               "chr",
                                               "pos"))
      
      if (!identical(colnames(geno), map$marker)) {
        stop("marker names in genotypic data and in map are not the same")
        
      }
    }
    
    
    # test that all genotypes present in the phenotypic data are also present in the genotypic data
    
    if (!(all(pheno$geno_ID %in% row.names(geno)))) {
      stop(
        "lines identified in the phenotypic data not identical to lines identified in the genotypic data"
      )
      
    }
    
    
    
    # test phenotypic data
    # test correct class for the different columns of the phenotype data
    
    
    if (!is.character(pheno[, 1])) {
      stop("the genotype names/IDs (first column of pheno) in pheno data must be character")
    }
    if (!is.numeric(pheno[, 2]) & !is.factor(pheno[, 2])) {
      stop("the year (second column of pheno) in pheno data must be numeric")
    }
    if (!is.character(pheno[, 3]) & !is.factor(pheno[, 3])) {
      stop("the location (third column of pheno) in pheno data must be character")
    }
    
    
    # Assign col.names pheno columns and transform year + location to factor
    
    pheno$year = as.factor(pheno$year)
    pheno$location = as.factor(pheno$location)
    
    # Give a numerical trait name if no name provided
    if (is.null(colnames(pheno)[4:ncol(pheno)])) {
      trait_names <- paste0("trait", 4:dim(pheno)[2])
      colnames(pheno) <- trait_names
      
    }
    
    
    
    # Check info_environments
    
    checkmate::assert_data_frame(info_environments,
                                 any.missing = F,
                                 min.cols = 4)
    class(info_environments) <- "data.frame"
    checkmate::assert_names(
      names(info_environments),
      must.include = c("year",
                       "location",
                       "longitude",
                       "latitude")
    )
    
    
    
    if (compute_climatic_ECs &
        is.null(info_environments$harvest.date)) {
      stop("Computation of ECs is required but no date for the harvest date.")
    }
    if (compute_climatic_ECs &
        is.null(info_environments$planting.date)) {
      stop("Computation of ECs is required but no date for the planting date.")
    }
    
    if (compute_climatic_ECs  &
        !inherits(info_environments$harvest.date, "Date")) {
      stop("planting date in info_environments as Date (format y-m-d).")
    }
    if (compute_climatic_ECs &
        !inherits(info_environments$planting.date, "Date")) {
      stop("harvest date in info_environments as Date (format y-m-d).")
    }
    
    
    # Create unique ID environment based on the location x year combination
    pheno$IDenv <- paste0(pheno$location, "_", pheno$year)
    info_environments$IDenv <-
      paste0(info_environments$location, "_", info_environments$year)
    
    
    # if geographical coordinates data.frame provided, test that all locations in the pheno data are present in the info_environments data.frame
    # test that longitude and latitude numerically provided
    #
    if (!is.data.frame(info_environments)) {
      stop("info_environments is not a data.frame")
    }
    
    if (!all(unique(pheno$IDenv) %in% unique(info_environments$IDenv))) {
      stop(
        "Environments identified in the info_environments data.frame are not identical to the locations present in the phenotypic data."
      )
      
    }
    if (!is.character(info_environments$location)) {
      stop("location is not character in info_environments")
    }
    
    if (!is.numeric(info_environments$year)) {
      stop("year is not numeric in info_environments")
    }
    
    if (!is.numeric(info_environments$longitude)) {
      stop("longitude is not numeric in info_environments")
    }
    
    if (!is.numeric(info_environments$latitude)) {
      stop("latitude is not numeric in info_environments")
    }
    
    
    
    
    # if marker data.frame provided, test marker names + chromosome info + positions provided
    
    if (!is.null(map)) {
      if (!is.character(map$marker)) {
        stop("the marker name (first column in map) must be character")
        
      }
      
      if (!is.numeric(map$chr)) {
        stop("the chromosome number (second column in map) must be numeric")
        
      }
      
      
      if (!is.numeric(map$pos)) {
        stop("the genetic position (third column in map) must be numeric")
        
      }
      
      
      
    } else {
      cat("No map provided.\n")
    }
    
    # test environmental data
    
    if (!is.null(climate_variables)) {
      if (nrow(climate_variables) != length(unique(pheno$IDenv))) {
        stop(
          "The number of observations in the climate_variables dataset does not match the number of Year x Location combinations from the pheno file."
        )
      }
      
      checkmate::assert_names(names(climate_variables),
                              must.include = c("year",
                                               "location"))
      
      
      if (!is.numeric(climate_variables[, 1])) {
        stop("The first column of climate_variables dataset should contain the year as numeric.")
      }
      if (!is.character(climate_variables[, 2])) {
        stop(
          "The second column of climate_variables dataset should contain the location as character."
        )
      }
      if (!all(
        vapply(
          climate_variables[, 3:ncol(climate_variables)],
          FUN = function(col) {
            is.numeric(col)
          },
          FUN.VALUE = logical(1),
          USE.NAMES = FALSE
        )
      )) {
        stop(
          "Col3+ of climate_variables dataset should contain the environmental variable as numeric."
        )
      }
      
      # Assign col.names of climate_variables
      
      
      climate_variables$IDenv <-
        paste0(climate_variables$location, "_", climate_variables$year)
      
      
      
    } else{
      cat("No climate covariates provided by the user.\n")
    }
    
    if (!is.null(soil_variables)) {
      if (nrow(soil_variables) != length(unique(pheno$IDenv))) {
        stop(
          "The number of observations in the soil variables dataset does not",
          "match the number of Year x Location combinations from the pheno",
          "file."
        )
      }
  
      
      if (!is.numeric(soil_variables[, 1])) {
        stop("The first column of soil_variables dataset should contain the year as numeric.")
      }
      if (!is.character(soil_variables[, 2])) {
        stop(
          "The second column of soil_variables dataset should contain the location as character."
        )
      }
      
      
      # Assign col.names of soil_variables
      
      checkmate::assert_names(names(soil_variables),
                              must.include = c("year",
                                               "location"))
      
      
      soil_variables$IDenv <-
        paste0(soil_variables$location, "_", soil_variables$year)
      
      
      
    } else{
      cat("No soil covariates provided by the user.\n")
      if (get_public_soil_data){
      cat("The package will retrieve soil data from the SoilGrids (ISRIC)",
          "Database.")
      soil_variables_list <- 
        lapply(
          unique(info_environments$IDenv),
          FUN = function(x, ...) {
            get_soil_per_env(
              environment = x,
              info_environments = info_environments,
              ...
            )
          }
        )
      
      soil_variables <- do.call("rbind", soil_variables_list)
      }
      
    }
    
    
    
    if (!compute_climatic_ECs & is.null(climate_variables)) {
      cat(
        paste(
          "No climate covariates will be computed nor used using",
          "the package. To allow calculation of ECs, please use the",
          "argument compute_climatic_ECs = T.\n"
        )
      )
    }
    
    if (compute_climatic_ECs &
        !is.null(climate_variables)) {
      stop(
        "Either climate variables (= environmental predictors) should be directly given, OR environmental predictors should be computed by the package. Raw weather daily data can be provided, according to the documentation."
      )
    }
    
    if (compute_climatic_ECs) {
      cat("Computation of environmental covariates starts.\n")
      
      merged_ECs <- get_ECs(
        info_environments = info_environments,
        raw_weather_data = raw_weather_data,
        path_data = path_to_save,
        ...
      )
      
      
      
      
      climate_variables <- merged_ECs$ECs
      climate_data_retrieved <- merged_ECs$climate_data_retrieved
      
      ECs_computed <- TRUE
      cat("Computation of environmental covariates is done.\n")
    }
    else{
      ECs_computed <- FALSE
      climate_data_retrieved <- FALSE
    }
    
    
    ### CLUSTERING OF ENVIRONMENTAL INFORMATION ###
    if (!is.null(path_to_save)) {
      if (!is.null(soil_variables) | !is.null(climate_variables)) {
        clustering_env_data(
          weather_ECs = climate_variables,
          soil_ECs = soil_variables,
          path_plots = path_to_save
        )
      }
    }
    
    ### MERGE climate_variables and soil_variables datasets
    if (!is.null(soil_variables) & !is.null(climate_variables)) {
      env_data <-
        merge(
          soil_variables %>% dplyr::select(-year,-location),
          climate_variables,
          by = c("IDenv")
        )
      list_climatic_predictors <-
        colnames(climate_variables %>% dplyr::select(-IDenv,-year,-location))
      list_soil_predictors <-
        colnames(soil_variables %>% dplyr::select(-IDenv,-year,-location))
    } else if (is.null(soil_variables) &
               !is.null(climate_variables)) {
      env_data <- climate_variables
      list_climatic_predictors <-
        colnames(climate_variables %>% dplyr::select(-IDenv,-year,-location))
      list_soil_predictors <- NULL
    } else if (!is.null(soil_variables) &
               is.null(climate_variables)) {
      env_data <- soil_variables
      list_climatic_predictors <- NULL
      list_soil_predictors <-
        colnames(soil_variables %>% dplyr::select(-IDenv,-year,-location))
      
    } else{
      env_data <- NULL
      list_climatic_predictors <- NULL
      list_soil_predictors <- NULL
    }
    
    
    
    
    METData <- structure(
      list(
        "geno" = geno,
        "map" = map,
        "pheno" = pheno,
        "ECs_computed" = ECs_computed,
        "climate_data_retrieved" = climate_data_retrieved,
        "env_data" = env_data,
        "info_environments" = info_environments,
        "list_climatic_predictors" = list_climatic_predictors,
        "list_soil_predictors" = list_soil_predictors
      ),
      class = c("METData", "list")
    )
    
    
    
    return(METData)
    
    
  }

#' @rdname create_METData
#' @aliases create_METData
#' @export
create_METData <- function(geno = NULL,
                           pheno = NULL,
                           info_environments = NULL,
                           map = NULL,
                           climate_variables = NULL,
                           compute_climatic_ECs = FALSE,
                           soil_variables = NULL,
                           raw_weather_data = NULL,
                           path_to_save = NULL,
                           ...) {
  validate_create_METData(
    new_create_METData(
      geno = geno,
      pheno = pheno,
      info_environments = info_environments,
      map = map,
      climate_variables = climate_variables,
      soil_variables = soil_variables,
      compute_climatic_ECs = compute_climatic_ECs,
      raw_weather_data = raw_weather_data,
      path_to_save = path_to_save,
      ...
    )
  )
}

#' @rdname create_METData
#' @aliases create_METData
#' @export
validate_create_METData <- function(x,
                                    ...) {
  checkmate::assert_class(x, "METData")
  
  checkmate::assert_names(
    names(x),
    must.include = c(
      "geno",
      "map",
      "pheno",
      "ECs_computed",
      "climate_data_retrieved",
      "env_data",
      "info_environments",
      "list_climatic_predictors",
      "list_soil_predictors"
      
    )
  )
  
  checkmate::assert_class(x[["geno"]], "matrix")
  checkmate::assertFALSE(checkmate::anyMissing(x[["geno"]]))
  
  checkmate::assert_data_frame(x[["map"]], null.ok = TRUE)
  
  checkmate::assert_class(x[["pheno"]], "data.frame")
  
  checkmate::assert_class(x[["env_data"]], "data.frame", null.ok = TRUE)
  
  checkmate::assert_class(x[["info_environments"]], "data.frame")
  checkmate::assertFALSE(checkmate::anyMissing(x[["info_environments"]]))
  
  
  
  checkmate::assert_class(x[["ECs_computed"]], "logical")
  checkmate::assertFALSE(checkmate::anyMissing(x[["ECs_computed"]]))
  
  checkmate::assert_class(x[["climate_data_retrieved"]], "logical")
  checkmate::assertFALSE(checkmate::anyMissing(x[["climate_data_retrieved"]]))
  
  
  checkmate::assert_character(x[["list_climatic_predictors"]], null.ok = TRUE)
  
  checkmate::assert_character(x[["list_soil_predictors"]], null.ok = TRUE)
  
  
  
  return(x)
}
