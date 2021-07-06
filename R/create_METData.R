#' Create a multi-environment trials data object
#'
#' @description
#' This function combines all types of data sources (genotypic, phenotypic,
#' environmental data, information about the environments...) in a single data
#' object of class \code{METData}.
#'
#' @param geno \code{numeric} genotype values stored in a \code{matrix} or
#'   \code{data.frame} which contains the geno_ID as row.names and markers as
#'   columns.
#'
#' @param pheno \code{data.frame} object with at least 4 columns.
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'   }
#'   From the fourth column on: each column \code{numeric} contains phenotypic
#'   values for a phenotypic trait observed in a combination Year x Location.
#'   Names of the traits can be provided as column names.
#'   \cr
#'   * \strong{The geno_ID must be the same as the row.names in the geno object.
#'   }
#'   \cr
#'   * \strong{For genotypes to be predicted for real case scenarios (e.g. no
#'   phenotypic data available), fill with NA in the column containing the
#'   phenotypic values of the trait.}
#'   \cr
#'
#' @param info_environments \code{data.frame} object with at least the 4 first
#'   columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'     \item planting.date: (optional) \code{Date} YYYY-MM-DD
#'     \item harvest.date: (optional) \code{Date} YYYY-MM-DD \cr
#'   }
#'   * \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location evaluated across four years, 4
#'   rows should be present.}
#'   * \strong{The fifth and sixth columns (planting.date and harvest.date) are
#'   required only if the user wants to download weather data with the
#'   package (setting argument `compute_ECs = T` in [create_METData()] and using
#'   subsequently the function [get_ECs()]).}
#'
#' @param map \code{data.frame} object with 3 columns.
#'   \enumerate{
#'   \item marker_name \code{character} with marker names
#'   \item chr \code{numeric} with chromosome number
#'   \item pos \code{numeric} with marker position.
#'   }
#'   \strong{Map object not absolutely required}.
#'
#' @param env_data \code{data.frame} can be let as NULL by user, if no
#'   environment data provided as input. Otherwise, a \code{data.frame} should
#'   be provided:
#'   Two types of \code{data.frame} can be provided:
#'   * **Type 1**: if `unique_EC_by_geno = FALSE`. Each environmental covariate
#'   characterizes all phenotypes within a complete environment (e.g. ECs are
#'   not computed individually for each variety within an environment).
#'   \strong{The data.frame should contain as many rows as
#'   the info_environments \code{data.frame}.} \cr
#'   Columns should contain:
#'   \enumerate{
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'   }
#'   Columns 3 and + should be numeric and contain the environmental covariates.
#'   \cr
#'
#'   OR \cr
#'   * **Type 2**: if `unique_EC_by_geno = TRUE`. Each environmental covariate
#'   is computed specifically for an environment AND for a genotype (e.g. ECs
#'   are computed individually for each variety within an environment).
#'   \strong{The data.frame should contain as many rows as the phenotypic
#'   dataset \code{data.frame}.} \cr
#'   Columns should contain:
#'   \enumerate{
#'     \item geno_ID \code{character} with genotype ID
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'   }
#'   Columns 4 and + should be numeric and contain the environmental covariates.
#'   \cr
#'
#'   * \strong{Providing env_data  and still setting `compute_ECs = T` is
#'   possible. For instance, the user can have some information regarding the
#'   soil composition (\% % clay, sand, silt, organic matter content).
#'   A disease status can also be encoded as categorical variable if it affected
#'   some environments. In addition to these covariates, weather-based
#'   covariates will be computed if `compute_ECs = T` with the package if they
#'   were not available to provide as input by the user.}
#'
#' @param unique_EC_by_geno \code{logical} indicates if the environmental
#'   covariates contained in env_data are also genotype-specific (dependent on
#'   the phenology, for instance) or unique for a complete environment. Default
#'   is `FALSE`.
#'
#' @param compute_ECs \code{logical} indicates if environmental covariates
#'   should be computed in further steps using the function [get_ECs()]. Default
#'   is `FALSE`. \cr
#'   \strong{Set compute_ECs = `TRUE` if user wants to use weather data
#'   acquired with this package.}
#'
#' @param filtering_markers \code{logical} indicator if a single-environment
#'   QTL detection step (performed with GWAS or ElasticNet using the function
#'   [select_markers()]) should be conducted, to identify a subset of markers
#'   representing potential environment-specific QTL effects, for which
#'   interactions with environmental covariates can be further investigated.
#'   Default is `TRUE`.
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
#'   covariates should be later computed.
#'
#' * **env_data**: \code{data.frame} with the environmental covariates per
#'   environment (and if genotype-specific, per genotype).
#'
#' * **info_environments**: \code{data.frame} contains basic information on
#'   each environment.
#'
#' * **unique_EC_by_geno**: \code{logical} to indicate if the EC is genotype-
#'   specific.
#'
#' * **filtering_markers**: \code{logical} indicates if a filtering marker
#'   step should be applied in further steps
#'
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' @examples
#'
#' data(geno_G2F)
#' data(pheno_G2F)
#' data(map_G2F)
#' data(info_environments_G2F)
#' METdata_G2F <- create_METData(geno=geno_G2F,pheno=pheno_G2F,map=map_G2F,env_data = NULL,compute_ECs = FALSE,info_environments = info_environments_G2F)
#'
#' data(geno_indica)
#' data(map_indica)
#' data(pheno_indica)
#' data(info_environments_indica)
#' data(env_data_indica)
#' METdata_indica <- create_METData(geno=geno_indica,pheno=pheno_indica,env_data = env_data_indica,unique_EC_by_geno = FALSE,compute_ECs = FALSE,info_environments = info_environments_indica,map = map_indica)
#'
#' data(geno_japonica)
#' data(map_japonica)
#' data(pheno_japonica)
#' pheno_japonica1<-pheno_japonica[-which(pheno_japonica$year==2013),]
#' data(info_environments_japonica)
#' info_environments_japonica1<-info_environments_japonica[-which(info_environments_japonica$year==2013),]
#' data(env_data_japonica)
#' env_data_japonica1<-env_data_japonica[-which(env_data_japonica$year==2013),]
#' METdata_japonica1 <- create_METData(geno=geno_japonica,pheno=pheno_japonica1,env_data = env_data_japonica1,unique_EC_by_geno = FALSE,compute_ECs = FALSE,info_environments = info_environments_japonica1,map = map_japonica)
#'




create_METData <-
  function(geno = NULL,
           pheno = NULL,
           info_environments = NULL,
           map = NULL,
           env_data = NULL,
           unique_EC_by_geno = FALSE,
           compute_ECs = FALSE,
           filtering_markers = TRUE,
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
    if (!is.numeric(pheno[, 2])) {
      stop("the year (second column of pheno) in pheno data must be numeric")
    }
    if (!is.character(pheno[, 3])) {
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
    
    # Assign col.names of env_data
    
    if (!is.null(env_data)) {
      if (unique_EC_by_geno == TRUE) {
        colnames(env_data)[1:3] <- c('geno_ID', 'year', 'location')
      }
      if (unique_EC_by_geno == FALSE) {
        colnames(env_data)[1:2] <- c('year', 'location')
      }
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
    if (compute_ECs &
        is.null(info_environments$harvest.date)) {
      stop('Computation of ECs is required but no date for the harvest date.')
    }
    if (compute_ECs &
        is.null(info_environments$planting.date)) {
      stop('Computation of ECs is required but no date for the planting date.')
    }
    
    if (compute_ECs &
        !inherits(info_environments_G2F$harvest.date, 'Date')) {
      stop('planting date in info_environments as Date (format y-m-d).')
    }
    if (compute_ECs &
        !inherits(info_environments_G2F$planting.date, 'Date')) {
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
      if (unique_EC_by_geno == FALSE &
          nrow(env_data) != length(unique(pheno$IDenv))) {
        stop(
          'The number of observations in the environmental data does not match the number of Year x Location combinations from the pheno file.'
        )
      }
      
      
      if (unique_EC_by_geno == TRUE &
          nrow(env_data) != nrow(pheno)) {
        stop(
          'The number of observations in the environmental data does not match the number of observations in the pheno file (Year x Location x Genotype).'
        )
      }
      
      
      
      ## Test specific format for environmental data: if genotype based environmental covariates (taking phenology into account to compute ECs)
      
      if (unique_EC_by_geno == TRUE) {
        if (!is.character(env_data[, 1])) {
          stop(
            'The first column of environmental data should contain the genotype names/IDs as character.'
          )
        }
        if (!is.numeric(env_data[, 2])) {
          stop('The second column of environmental data should contain the year as numeric.')
        }
        if (!is.character(env_data[, 3])) {
          stop(
            'The third column of environmental data should contain the location as character.'
          )
        }
        if (!all(vapply(
          env_data[, 4:ncol(env_data)],
          FUN = function(col) {
            is.numeric(col)
          },
          FUN.VALUE = logical(1),
          USE.NAMES = FALSE
        ))) {
          stop(
            'Col4+ of environmental data should contain the environmental variable as numeric.'
          )
        }
      }
      if (unique_EC_by_geno == FALSE) {
        if (!is.numeric(env_data[, 1])) {
          stop('The first column of environmental data should contain the year as numeric.')
        }
        if (!is.character(env_data[, 2])) {
          stop(
            'The second column of environmental data should contain the location as character.'
          )
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
      }
      
      env_data$IDenv <-
        paste0(env_data$location, '_', env_data$year)
      
      
      
    } else{
      cat('No environmental covariates provided by the user.\n')
    }
    
   
    
    if (!compute_ECs & is.null(env_data)) {
      cat(
        paste(
          'No environmental covariates will be computed nor used using',
          'the package. To allow calculation of ECs, please use the',
          'argument compute_ECs = T.\n'
        )
      )
    }
    
    if (compute_ECs) {
      cat('Computation of environmental covariates starts.\n')
      env_data <- get_ECs(info_environments = info_environments,
                          ...)
      
      # Add ECs to the table env_data, if this table already contains 
      # environmental covariates.
      
      if (!is.null(env_data) & !unique_EC_by_geno) {
        env_data <-
          merge(merged_ECs, env_data %>%select(-location,-year), by = 'IDenv')
      }
      
      
      
      else{
        env_data <- merged_ECs
      }
      
      ECs_computed <- TRUE
      cat('Computation of environmental covariates is done.\n')
      METData <- structure( list(
        'geno' = geno,
        'map_markers' = map,
        'pheno' = pheno,
        'compute_ECs' = compute_ECs,
        'env_data' = env_data,
        'info_environments' = info_environments,
        'unique_EC_by_geno' = unique_EC_by_geno,
        'filtering_markers' = filtering_markers,
        'ECs_computed' = ECs_computed
      ), class='METData')
    }
    else{
    
    
    METData <- structure( list(
      'geno' = geno,
      'map_markers' = map,
      'pheno' = pheno,
      'compute_ECs' = compute_ECs,
      'env_data' = env_data,
      'info_environments' = info_environments,
      'unique_EC_by_geno' = unique_EC_by_geno,
      'filtering_markers' = filtering_markers
    ), class='METData')}
    
    return(METData)
    
    
  }
