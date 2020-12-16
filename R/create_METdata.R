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
#' @param
#'
#' @param info_environments \code{data.frame} with at least 4 columns
#'   First column \code{numeric} with the year
#'   Second column \code{Character} with the location
#'   Third column \code{numeric} with the longitude
#'   Fourth column \code{numeric} with the latitude
#'
#' @backref
#'
#'
#'
#'
#'
#'


# repeated location lines if same year
#  columns in info_environments_
'year', 'Location', 'Longitude', 'Latitude', 'planting.date', 'Harvest.Date'


## missing phenotypes for genotypes to predict should be indicated by NA


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



    # Match the coordinates information with the environmental data






    METpred_data <- list(
      'geno' = geno,
      'pheno' = pheno,
      'compute_ECs' = compute_ECs,
      'env_data' = env_data,
      'map_markers' = map,
      'filtering_markers' = filtering_markers
    )
    return(METpred_data)


  }


object <-
  create_METdata(
    geno = geno,
    pheno = pheno,
    env_data = env_data,
    unique_EC_by_geno = TRUE,
    compute_ECs = FALSE,
    info_environments = info_environments
  )
