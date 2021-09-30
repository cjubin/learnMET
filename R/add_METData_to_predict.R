#' Create a multi-environment trials data object, which will be used for 
#' prediction of unobserved data.
#'
#' @description
#' This function combines all types of predictor sources (genotypic,
#' information about the environments to be predicted, environmental data if 
#' available...) and the data.frame describing the phenotypic observations to 
#' be predicted. 
#' 
#' @name add_METData_to_predict
#' 
#' @param METData_training \code{METData} object, which will be used as training
#'   set. Normally, it should be the same \code{METData} object as the one 
#'   created with [create_METData()] and evaluated with [predict_trait_MET_cv()].
#'   
#' @param geno_new \code{data.frame} object which should contain genotypic data
#'   for new candidates, if these are not included in `METData_training$geno`
#'   which will be used as training set. If already included, no need to include
#'   geno_new. Default is `NULL`.
#'
#' @param pheno_new \code{data.frame} object with at least 3 columns.
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'   }
#'   Basically, this `data.frame` indicates which genotypes should be
#'   predicted in which environments. No columns corresponding to phentoypic 
#'   traits is expected.
#'
#' @param info_environments_to_predict \code{data.frame} object with at least 
#'   the 4 following columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'  }
#'  The two next columns are required only if weather data should be 
#'  retrieved from NASA POWER data using the argument `compute_climatic_EC` set 
#'  to TRUE.
#'  \enumerate{
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD \cr
#'   }
#'   * \strong{The data.frame should contain as many rows as Year x Location
#'   combinations which will be used in pheno_new.}
#' 
#' @param env_data_manual \code{data.frame} can be let as NULL by user, if no
#'   environmental covariables provided as input. Otherwise, a \code{data.frame}
#'   should be provided.
#'   \strong{The data.frame should contain as many rows as the phenotypic
#'   dataset \code{data.frame}.} \cr
#'   Columns should be:
#'   \enumerate{
#'     \item geno_ID \code{character} with genotype ID
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'   }
#'   Columns 4 and + should be numeric and contain the environmental covariates.
#'   \cr
#'
#'   * \strong{Providing env_data  and setting `compute_climatic_ECs` to `TRUE`
#'   is possible. For instance, the user can have some information regarding the
#'   soil composition (\% % clay, sand, silt, organic matter content).
#'   A disease status can also be encoded as categorical variable if it affects
#'   some environments. In addition to these type of covariates, weather-based
#'   covariates will be computed if `compute_climatic_ECs` is set to `TRUE`.}
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
#' @param compute_climatic_ECs \code{logical} indicates if climatic covariates 
#'   should be computed with the function. Default is `FALSE`. \cr
#'   \strong{Set compute_climatic_ECs = `TRUE` if user wants to use weather data
#'   from NASA POWER data. For instance, if no weather-based covariables
#'   can be provided or if raw weather data are only available for some 
#'   environments but not for others.}
#' 
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
#' * **compute_climatic_ECs**: \code{logical} indicates if environmental
#'   covariates were required to be retrieved via the package by the user.
#'
#' * **env_data**: \code{data.frame} with the environmental covariates per
#'   environment.
#'
#' * **info_environments**: \code{data.frame} contains basic information on
#'   each environment.
#'   
#' * **ECs_computed**: OPTIONAL \code{logical} subelement added in the output
#'   only if the function [get_ECs()] was correctly run within the pipeline.
#'
#' @author Cathy C. Westhues \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' @examples
#'
new_add_METData_to_predict <-
  function(METData_training,
           geno_new = NULL,
           pheno_new = NULL,
           info_environments_to_predict = NULL,
           env_data_manual = NULL,
           raw_weather_data = NULL,
           compute_climatic_ECs = FALSE,
           ...) {
    
    ################################
    ## Check object pheno_new ##
    ################################
    
    
    if (is.null(pheno_new)) {
      stop("The table of phenotypic observations to predict must be provided.")
    }
    
    
    ## Check the pheno_new input data ##
    
    if (!is.data.frame(pheno_new)) {
      stop('pheno_new should be given as a data.frame')
    }
    if (ncol(pheno_new) < 3) {
      stop(
        'pheno_new_data should contain 3 columns: genotype lines (col1), year (col2), location (col3).'
      )
    }
    
    if (!is.character(pheno_new[, 1])) {
      stop("the genotype names/IDs (first column of pheno_new) must be character")
    }
    if (!is.numeric(pheno_new[, 2]) & !is.factor(pheno_new[, 2])) {
      stop("the year (second column of pheno_new) must be numeric or factor.")
    }
    if (!is.character(pheno_new[, 3]) &
        !is.factor(pheno_new[, 3])) {
      stop("the location (third column of pheno_new) must be character or factor.")
    }
    
    
    # Assign col.names to pheno_new columns and transform year + location to
    # factor
    
    colnames(pheno_new)[1:3] <- c('geno_ID', 'year', 'location')
    pheno_new$year = as.factor(pheno_new$year)
    pheno_new$location = as.factor(pheno_new$location)
    
    
    # Create unique ID environment based on the location x year combination
    pheno_new$IDenv <-
      paste0(pheno_new$location, '_', pheno_new$year)
    
    
    
    
    
    ################################
    ## PREPARE ENVIRONMENTAL DATA ##
    ################################
    
    
    if (compute_climatic_ECs &
        is.null(info_environments_to_predict)) {
      stop(
        paste(
          "Basic information is needed about the enviroments for which",
          "environmental data should be retrieved."
        )
      )
    }
    
    if (compute_climatic_ECs &
        !is.null(info_environments_to_predict)) {
      if (ncol(info_environments_to_predict) < 6) {
        stop(
          'info_environments_to_predict should contain at least 6 columns: year, location, longitude, latitude, planting.date and harvest.date (in this order).'
        )
      }
      colnames(info_environments_to_predict)[1:6] <-
        c('year',
          'location',
          'longitude',
          'latitude',
          'planting.date',
          'harvest.date')
      if (compute_climatic_ECs &
          is.null(info_environments_to_predict$harvest.date)) {
        stop('Computation of ECs is required but no date for the harvest date.')
      }
      if (compute_climatic_ECs &
          is.null(info_environments_to_predict$planting.date)) {
        stop('Computation of ECs is required but no date for the planting date.')
      }
      
      if (compute_climatic_ECs &
          !inherits(info_environments_to_predict$harvest.date, 'Date')) {
        stop('planting date in info_environments_to_predict as Date (format y-m-d).')
      }
      if (compute_climatic_ECs &
          !inherits(info_environments_to_predict$planting.date, 'Date')) {
        stop('harvest date in info_environments_to_predict as Date (format y-m-d).')
      }
      
      info_environments_to_predict$IDenv <-
        paste0(
          info_environments_to_predict$location,
          '_',
          info_environments_to_predict$year
        )
      
    }
    
    if (!all(unique(pheno_new$IDenv) %in% unique(info_environments_to_predict$IDenv))) {
      stop(
        "Environments identified in the info_environments_to_predict data.frame are not identical to the locations present in the phenotypic data."
      )
      
    }
    if (!is.character(info_environments_to_predict$location)) {
      stop("location is not character in info_environments")
    }
    
    if (!is.numeric(info_environments_to_predict$year)) {
      stop("year is not numeric in info_environments")
    }
    
    if (!is.numeric(info_environments_to_predict$longitude)) {
      stop("longitude is not numeric in info_environments")
    }
    
    if (!is.numeric(info_environments_to_predict$latitude)) {
      stop("latitude is not numeric in info_environments")
    }
    
    
    
    if (compute_climatic_ECs | !is.null(raw_weather_data)) {
      cat('Computation of environmental covariates starts for untested environments.\n')
      merged_ECs <- get_ECs(info_environments = info_environments,
                            raw_weather_data = raw_weather_data,
                            ...)
      
      
      # Add ECs to the table env_data, if this table already contains
      # environmental covariates.
      
      if (!is.null(env_data_manual)) {
        env_data_new <-
          merge(merged_ECs,
                env_data_manual %>% select(-location,-year),
                by = 'IDenv')
      }
      
      
      
      else{
        env_data_new <- merged_ECs
      }
      
      ECs_computed <- TRUE
      
      
      cat('Computation of environmental covariates for environments to predict is done.\n')
      
    } else {
      ECs_computed <- FALSE
    }
    
    
    ################################
    ## PREPARE GENOMIC DATA ##
    ################################
    
    
    ## Test the geno_new data input ##
    
    if (is.null(geno_new) &
        all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
      geno_new <- METData_training$geno
      map_new <- METData_training$map
      
    } else if (is.null(geno_new) &
               !all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
      stop(
        'Some genotype names are not present in the row.names of METData_training$geno, and no additional genotype data were provided.
           Please add geno data for the lines to be predicted which are not in row.names(METData_training$geno)'
      )
    } else {
      # test format of the genotypic data
      
      if (!is.matrix(geno_new) & !is.data.frame(geno_new)) {
        stop("genotypic data not provided as a matrix or data.frame")
      } else{
        geno_new <- as.data.frame(geno_new)
      }
      
      
      if (!all(apply(geno_new, 2, is.numeric))) {
        stop("genotypic data not provided as numeric")
      }
      
      
      
      common_cols <-
        intersect(colnames(METData_training$geno), colnames(geno_new))
      geno_new <- rbind(METData_training$geno[, common_cols],
                        geno_new[, common_cols])
      geno_new$geno_ID <- row.names(geno_new)
      geno_new <- unique(geno_new) %>% dplyr::select(-geno_ID)
      
      map_new <-
        METData_training$map[which(METData_training$map$marker_name %in%
                                             common_cols),]
      
    }
    
    
    
    ################################
    ## CHECK/PREPARE ENV DATA ##
    ################################
    
    if (!is.null(env_data_manual)) {
      if (nrow(env_data_manual) != length(unique(pheno_new$IDenv))) {
        stop(
          'The number of observations in the environmental data does not match the number of Year x Location combinations from the pheno file.'
        )
      }
      
      
      ## Test specific format for environmental data: if genotype based environmental covariates (taking phenology into account to compute ECs)
      
      
      
      if (!is.numeric(env_data_manual[, 1]) &
          !is.factor(env_data_manual[, 2])) {
        stop('The first column of environmental data should contain the year as numeric.')
      }
      if (!is.character(env_data_manual[, 2]) &
          !is.factor(env_data_manual[, 3])) {
        stop('The second column of environmental data should contain the location as character.')
      }
      if (!all(vapply(
        env_data_manual[, 3:ncol(env_data_manual)],
        FUN = function(col) {
          is.numeric(col)
        },
        FUN.VALUE = logical(1),
        USE.NAMES = FALSE
      ))) {
        stop(
          'Col3+ of environmental data should contain the environmental variable as numeric.'
        )
      }
      
      
      env_data_manual$IDenv <-
        paste0(env_data_manual$location, '_', env_data_manual$year)
      
      
      env_data_new <- merge(
        env_data_new,
        info_environments_to_predict %>% dplyr::select(-year,-location,-longitude,-latitude),
        by = 'IDenv',
        all.x = T
      )
      
    }
    
    
  
      METData <- structure(
        list(
          'geno' = geno_new,
          'map' = map_new,
          'pheno' = pheno_new,
          'compute_climatic_ECs' = compute_climatic_ECs,
          'ECs_computed' = ECs_computed,
          'env_data' = env_data_new,
          'info_environments' = info_environments_to_predict
        ),
        class = c('METData','list')
      )
      
      return(METData)
    
    
    
  }






#' @rdname add_METData_to_predict
#' @aliases add_METData_to_predict
#' @export
add_METData_to_predict <- function(METData_training,
                            pheno_new = NULL,
                            geno_new = NULL,
                            compute_climatic_ECs = FALSE,
                            info_environments_to_predict = NULL,
                            env_data_manual = NULL,
                            raw_weather_data = NULL,
                            ...) {
  
  validate_add_METData_to_predict(
    new_add_METData_to_predict(
      METData_training = METData_training,
      pheno_new = pheno_new,
      geno_new = geno_new,
      compute_climatic_ECs = compute_climatic_ECs,
      info_environments_to_predict = info_environments_to_predict,
      env_data_manual = env_data_manual,
      raw_weather_data = raw_weather_data,
      ...
    )
  )
}


#' @rdname add_METData_to_predict
#' @aliases add_METData_to_predict
#' @export

validate_add_METData_to_predict <- function(x,...) {
  
  checkmate::assert_class(x, 'METData')
  
  checkmate::assert_names(names(x), must.include = c('geno','map','pheno','compute_climatic_ECs','ECs_computed','env_data','info_environments'))

  
  checkmate::assert_class(x[['geno']], 'data.frame')
  checkmate::assertFALSE(checkmate::anyMissing(x[['geno']]))
  
  checkmate::assert_data_frame(x[['map']],null.ok = TRUE)
 
  checkmate::assert_class(x[['pheno']], 'data.frame')
  
  checkmate::assert_class(x[['env_data']], 'data.frame')
  checkmate::assertFALSE(checkmate::anyMissing(x[['env_data']]))
  
  checkmate::assert_class(x[['info_environments']], 'data.frame')
  checkmate::assertFALSE(checkmate::anyMissing(x[['info_environments']]))
  
  checkmate::assert_class(x[['compute_climatic_ECs']], 'logical')
  checkmate::assertFALSE(checkmate::anyMissing(x[['compute_climatic_ECs']]))
  
  checkmate::assert_class(x[['ECs_computed']], 'logical')
  checkmate::assertFALSE(checkmate::anyMissing(x[['ECs_computed']]))
  
  
  
  
  return(x)
}
