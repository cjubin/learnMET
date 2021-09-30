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
#'   \item marker_name \code{character} with marker names
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
#'  to TRUE, or if raw weather data are provided.
#'  \enumerate{
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD \cr
#'   }
#'   * \strong{The data.frame should contain as many rows as Year x Location
#'   combinations which will be used in pheno_new.}
#'  
#' @param env_data \code{data.frame} can be let as NULL by user, if no
#'   environment data provided as input. Otherwise, a \code{data.frame} should
#'   be provided.
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
#'
#' @param compute_climatic_ECs \code{logical} indicates if climatic covariates 
#'   should be computed with the function. Default
#'   is `FALSE`. \cr
#'   \strong{Set compute_climatic_ECs = `TRUE` if user wants to use weather data
#'   from NASA POWER data. For instance, if no weather-based covariables
#'   can be provided or if raw weather data are only available for some 
#'   environments but not for others.}
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
#'   environment (and if genotype-specific, per genotype).
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
#' data(geno_G2F)
#' data(pheno_G2F)
#' data(map_G2F)
#' data(info_environments_G2F)
#' METdata_G2F <- create_METData(geno=geno_G2F,pheno=pheno_G2F,map=map_G2F,env_data = NULL,compute_climatic_ECs = FALSE,info_environments = info_environments_G2F)
#'
#' data(geno_indica)
#' data(map_indica)
#' data(pheno_indica)
#' data(info_environments_indica)
#' data(env_data_indica)
#' METdata_indica <- create_METData(geno=geno_indica,pheno=pheno_indica,env_data = env_data_indica,compute_climatic_ECs = FALSE,info_environments = info_environments_indica,map = map_indica)
#'
#' data(geno_japonica)
#' data(map_japonica)
#' data(pheno_japonica)
#' pheno_japonica1<-pheno_japonica[-which(pheno_japonica$year==2013),]
#' data(info_environments_japonica)
#' info_environments_japonica1<-info_environments_japonica[-which(info_environments_japonica$year==2013),]
#' data(env_data_japonica)
#' env_data_japonica1<-env_data_japonica[-which(env_data_japonica$year==2013),]
#' METdata_japonica1 <- create_METData(geno=geno_japonica,pheno=pheno_japonica1,env_data = env_data_japonica1,compute_climatic_ECs = FALSE,info_environments = info_environments_japonica1,map = map_japonica)

new_create_METData <-
  function(geno = NULL,
           map = NULL,
           pheno = NULL,
           info_environments = NULL,
           env_data = NULL,
           raw_weather_data = NULL,
           compute_climatic_ECs = FALSE,
           path_to_save = NULL,
           ...) {
    # check if one object is missing
    
    if (is.null(geno)) {
      stop("genotypic data not provided")
    }
    
    if (is.null(pheno)) {
      stop("phenotypic data not provided")
    }
    
    if (is.null(info_environments)) {
      stop("info_environments not provided")
    }
    
    # test format of the genotypic data
    
    if (!is.matrix(geno) & !is.data.frame(geno)) {
      stop("genotypic data not provided as a matrix or data.frame")
    } else{
      geno <- as.data.frame(geno)
    }
    
    
    if (!all(apply(geno, 2, is.numeric))) {
      stop("genotypic data not provided as numeric")
    }
    
    # test that all genotypes present in the phenotypic data are also present in the genotypic data
    
    if (!(all(pheno[, 1] %in% row.names(geno)))) {
      stop(
        "lines identified in the phenotypic data not identical to lines identified in the genotypic data"
      )
      
    }
    
    # if marker matrix is given, test that the marker names are the same in the map and in the marker genotype matrices
    
    if (!is.null(map) & !identical(colnames(geno), map[, 1])) {
      stop("marker names in genotypic data and in map are not the same")
      
    }
    
    # test phenotypic data
    # test correct class for the different columns of the phenotype data
    if (!is.data.frame(pheno)) {
      
    }
    if (ncol(pheno) < 4) {
      stop(
        'MET pheno data should contain at least 4 columns: genotype lines (col1), year (col2), location (col3) and phenotypic values for at least one trait (from col4)'
      )
    }
    
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
    
    colnames(pheno)[1:3] <- c('geno_ID', 'year', 'location')
    pheno$year = as.factor(pheno$year)
    pheno$location = as.factor(pheno$location)
    
    # Give a numerical trait name if no name provided
    if (is.null(colnames(pheno)[4:ncol(pheno)])) {
      trait_names <- paste0('trait', 4:dim(pheno)[2])
      colnames(pheno) <- trait_names
      
    }
    
    
    
    # Assign colnames for info_environments columns
    if (ncol(info_environments) < 4) {
      stop(
        'info_environments should contain at least 4 columns: year, location, longitude, latitude.'
      )
    }
    colnames(info_environments)[1:4] <-
      c('year',
        'location',
        'longitude',
        'latitude')
    if ((compute_climatic_ECs|!is.null(raw_weather_data)) &
        is.null(info_environments$harvest.date)) {
      stop('Computation of ECs is required but no date for the harvest date.')
    }
    if ((compute_climatic_ECs | !is.null(raw_weather_data)) &
        is.null(info_environments$planting.date)) {
      stop('Computation of ECs is required but no date for the planting date.')
    }
    
    if ((compute_climatic_ECs | !is.null(raw_weather_data)) &
        !inherits(info_environments$harvest.date, 'Date')) {
      stop('planting date in info_environments as Date (format y-m-d).')
    }
    if ((compute_climatic_ECs | !is.null(raw_weather_data)) &
        !inherits(info_environments$planting.date, 'Date')) {
      stop('harvest date in info_environments as Date (format y-m-d).')
    }
    
    
    
    # Create unique ID environment based on the location x year combination
    pheno$IDenv <- paste0(pheno$location, '_', pheno$year)
    info_environments$IDenv <-
      paste0(info_environments$location, '_', info_environments$year)
    
    
    
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
      if (!is.character(map[, 1])) {
        stop("the marker name (first column in map) must be character")
        
      }
      
      if (!is.numeric(map[, 2])) {
        stop("the chromosome number (second column in map) must be numeric")
        
      }
      
      
      if (!is.numeric(map[, 3])) {
        stop("the genetic position (third column in map) must be numeric")
        
      }
      
      colnames(map) <- c('marker_name', 'chr', 'pos')
      
    } else {
      cat('No map provided.\n')
    }
    
    # test environmental data
    
    if (!is.null(env_data)) {
      if (nrow(env_data) != length(unique(pheno$IDenv))) {
        stop(
          'The number of observations in the environmental data does not match the number of Year x Location combinations from the pheno file.'
        )
      }
      
      
      ## Test specific format for environmental data: if genotype based environmental covariates (taking phenology into account to compute ECs)
      
      
      
      if (!is.numeric(env_data[, 1])) {
        stop('The first column of environmental data should contain the year as numeric.')
      }
      if (!is.character(env_data[, 2])) {
        stop('The second column of environmental data should contain the location as character.')
      }
      if (!all(vapply(
        env_data[, 3:ncol(env_data)],
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
      
      # Assign col.names of env_data
      
      colnames(env_data)[1:2] <- c('year', 'location')
      
      env_data$IDenv <-
        paste0(env_data$location, '_', env_data$year)
      
      
      
    } else{
      cat('No environmental covariates provided by the user.\n')
    }
    
    
    
    if (!compute_climatic_ECs & is.null(env_data)) {
      cat(
        paste(
          'No environmental covariates will be computed nor used using',
          'the package. To allow calculation of ECs, please use the',
          'argument compute_climatic_ECs = T.\n'
        )
      )
    }
    
    if (compute_climatic_ECs | !is.null(raw_weather_data)) {
      cat('Computation of environmental covariates starts.\n')
      merged_ECs <- get_ECs(info_environments = info_environments,
                            raw_weather_data = raw_weather_data,
                            path_data = path_to_save,
                            ...)
      
      # Add ECs to the table env_data, if this table already contains
      # environmental covariates.
      
      if (!is.null(env_data)) {
        env_data <-
          merge(merged_ECs, env_data %>% select(-location,-year), by = 'IDenv')
      }
      
      
      
      else{
        env_data <- merged_ECs
      }
      
      ECs_computed <- TRUE
      cat('Computation of environmental covariates is done.\n')}
      else{ECs_computed <- FALSE}
      
      
      
      
      METData <- structure(
        list(
          'geno' = geno,
          'map' = map,
          'pheno' = pheno,
          'compute_climatic_ECs' = compute_climatic_ECs,
          'ECs_computed' = ECs_computed,
          'env_data' = env_data,
          'info_environments' = info_environments,
          'ECs_computed' = ECs_computed
        ),
        class = c('METData','list')
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
                           env_data = NULL,
                           compute_climatic_ECs = FALSE,
                           raw_weather_data = NULL,
                           ...) {
  
  validate_create_METData(
    new_create_METData(
      geno = geno,
      pheno = pheno,
      info_environments = info_environments,
      map = map,
      env_data = env_data,
      compute_climatic_ECs = compute_climatic_ECs,
      raw_weather_data = raw_weather_data,
      ...
    )
  )
}

#' @rdname create_METData
#' @aliases create_METData
#' @export
validate_create_METData <- function(x,
                                    ...) {
  
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
