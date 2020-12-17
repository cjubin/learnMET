#' Create a multi-environment trials data object
#'
#'
#' @param geno \code{numeric} genotype values stored in a \code{matrix} which
#'   contains the geno_ID as row.names and markers as columns.
#'
#' @param pheno \code{data.frame} with at least 4 columns.
#'   First column "geno_ID"  \code{Character} contains the genotype identifiers.
#'   \strong{The geno_ID must be the same to the row.names in the geno object}
#'   \strong{For genotypes to be predicted (only geno data, no pheno values),
#'   fill with NA in pheno}
#'   Second column "year"  \code{numeric} contains the year for the observation
#'   Third column "location" \code{Character} contains the name of the location
#'   Fourth column and + \code{numeric} contain the phenotypic values for
#'   different traits. Names of the traits can be provided as col.names.
#'
#' @param info_environments \code{data.frame} with at least 4 columns.
#'   First column \code{numeric} with the year label
#'   Second column \code{Character} with the location
#'   Third column \code{numeric} with the longitude
#'   Fourth column \code{numeric} with the latitude
#'   Additional columns 'planting.date' and 'harvest.date' can be provided if
#'   these dates should be used to derive environmental covariates.
#'   \strong{The data.frame should contain as many rows as Year x Location
#'   combinations. Example: if only one location used in MET for several years,
#'   use still 1 row for each year.}
#'
#' @param map \code{data.frame} with 3 columns.
#'   First column \code{character} with marker names
#'   Second column \code{numeric} with chromosome number
#'   Third column \code{numeric} with marker position.
#'   \strong{not absolutely required}.
#'
#' @param env_data \code{data.frame} can be let as NULL by user, if no
#'   environment data already present.
#'   If env_data provided, and if unique_EC_by_geno==TRUE, the data.frame should
#'   contain as many rows as the pheno data.frame.
#'   First column \code{character} with genotype identifiers
#'   Second column \code{numeric} with the year label
#'   Third column \code{character} with the location character
#'   Columns 4+ should be numeric and contain the environmental covariates
#'   provided by the user.
#'   If unique_EC_by_geno==FALSE, the data.frame should contain as many rows as
#'   the info_environments data.frame:
#'   First column \code{numeric} with the year label
#'   Second column \code{character} with the location character
#'   Columns 3+ should be numeric and contain the environmental covariates
#'
#' @param unique_EC_by_geno \code{Logical} indicates if the environmental
#'  covariates contained in env_data are also genotype-specific (dependent on
#'  the phenology, for instance) or unique for a whole environment.
#'
#' @param compute_ECs \code{Logical} indicates if environmental covariates
#'   should be computed in further steps.
#'
#' @param filtering_markers \code{Logical} indicator if the number of markers
#' to use in the machine-learning based predictions should be pre-handled
#' in further steps (elasticnet, PCA...)
#'
#' @return
#'
#' a \code{list} of class \code{METData} which contains the following elements
#'
#' \item{geno}{\code{matrix} with genotype values of phenotyped individuals.}
#'
#' \item{map}{\code{data.frame} with genetic map.}
#'
#' \item{pheno}{\code{data.frame} with phenotypic trait values.}
#'
#' #' \item{compute_EC_by_geno}{\code{Logical} indicates if environmental
#' covariates should be later computed.}
#'
#' \item{env_data}{\code{data.frame} with the environmental covariates per
#' environment (and if genotype-specific, per genotype).}
#'
#' \item{info_environments}{\code{data.frame} contains basic information on
#' each environment.}
#'
#' \item{unique_EC_by_geno}{\code{Logical} to indicate if the EC is genotype-
#' specific.}
#'
#' \item{filtering_markers}{\code{Logical} indicates if a filtering marker step
#' should be applied in further steps}
#'
#' @examples
#'
#' data(geno_G2F)
#' data(pheno_G2F)
#' data(info_environments_G2F)
#' data(env_data_G2F)
#' METdata_G2F <- create_METdata(geno=geno_G2F,pheno=pheno_G2F,env_data = env_data_G2F,unique_EC_by_geno = T,compute_ECs = F,info_environments = info_environments_G2F)
#'
#' data(geno_indica)
#' data(map_indica)
#' data(pheno_indica)
#' data(info_environments_indica)
#' data(env_data_indica)
#' METdata_indica <- create_METdata(geno=geno_indica,pheno=pheno_indica,env_data = env_data_indica,unique_EC_by_geno = F,compute_ECs = F,info_environments = info_environments_indica,map = map_indica)
#'
#' data(geno_japonica)
#' data(map_japonica)
#' data(pheno_japonica)
#' data(info_environments_japonica)
#' data(env_data_japonica)
#' METdata_japonica <- create_METdata(geno=geno_japonica,pheno=pheno_japonica,env_data = env_data_japonica,unique_EC_by_geno = F,compute_ECs = F,info_environments = info_environments_japonica,map = map_japonica)
#'



create_METdata <-
  function(geno = NULL,
           pheno = NULL,
           map = NULL,
           env_data = NULL,
           unique_EC_by_geno = FALSE,
           compute_ECs = FALSE,
           info_environments = NULL,
           filtering_markers = TRUE) {
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

    if (!is.data.frame(geno)) {
      geno <- as.matrix(geno)
    }


    if (!is.matrix(geno)) {
      stop("genotypic data not provided as a matrix")
    }

    if (!is.numeric(geno)) {
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
        'MET pheno data should contain at least 4 columns: genotype lines (col1), year (col2), location (col3) and phenotype values (from col4)'
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

    # Assign col.names for info_environments columns
    if (ncol(info_environments)<4){stop('info-environments should contain at least 4 columns: year, location, longitude, latitude.')}
    colnames(info_environments)[1:4] <-
      c('year',
        'location',
        'longitude',
        'latitude')



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

    if (!identical(unique(as.character(as.vector(
      info_environments$IDenv
    ))), unique(as.character(as.vector(pheno$IDenv))))) {
      stop(
        "locations identified in the geographical coordinates data.frame are not identical to the locations present in the phenotypic data."
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

     colnames(map)<-c('marker_name','chr','pos')

    } else {
      cat('No map provided')
    }

    # test environmental data

    if (!is.null(env_data)) {
      if (unique_EC_by_geno == FALSE &&
          nrow(env_data) != length(unique(pheno$IDenv))) {
        stop(
          'The number of observations in the environmental data does not match the number of Year x Location combinations from the pheno file.'
        )
      }


      if (unique_EC_by_geno == TRUE &&
          nrow(env_data) != nrow(pheno)) {
        stop(
          'The number of observations in the environmental data does not match the number of observations in the pheno file (Year x Location x Genotype).'
        )
      }



      ## Test specific format for environmental data: if genotype based environmental covariates (taking phenology into account to compute ECs)

      if (unique_EC_by_geno == TRUE &&
          !is.character(env_data[, 1])) {
        stop(
          'The first column of environmental data should contain the genotype names/IDs as character.'
        )
      }
      if (unique_EC_by_geno == TRUE &&
          !is.numeric(env_data[, 2])) {
        stop('The second column of environmental data should contain the year as numeric.')
      }
      if (unique_EC_by_geno == TRUE &&
          !is.character(env_data[, 3])) {
        stop('The third column of environmental data should contain the location as character.')
      }
      if (unique_EC_by_geno == TRUE &&
          !all(vapply(
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

      if (unique_EC_by_geno == FALSE &&
          !is.numeric(env_data[, 1])) {
        stop('The first column of environmental data should contain the year as numeric.')
      }
      if (unique_EC_by_geno == FALSE &&
          !is.character(env_data[, 2])) {
        stop('The second column of environmental data should contain the location as character.')
      }
      if (unique_EC_by_geno == FALSE &&
          !all(vapply(
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

      env_data$IDenv <-
        paste0(env_data$location, '_', env_data$year)


      env_data <- merge(
        env_data,
        info_environments %>% select(-year, -location),
        by = 'IDenv',
        all.x = T
      )

    } else{
      cat('\nNo environmental covariates provided')
    }
    if (compute_ECs == TRUE) {
      cat('\nEnvironmental covariates should be determined.')
    }




    METData <- list(
      'geno' = geno,
      'map_markers' = map,
      'pheno' = pheno,
      'compute_ECs' = compute_ECs,
      'env_data' = env_data,
      'info_environments' = info_environments,
      'unique_EC_by_geno' = unique_EC_by_geno,
      'filtering_markers' = filtering_markers
    )

    class(METData)<- c("METData", "list")
    return(METData)


  }
