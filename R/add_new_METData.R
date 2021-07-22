#' Create a multi-environment trials data object
#'
#' @description
#' This function combines all types of data sources (genotypic, phenotypic,
#' environmental data, information about the environments...) in a single data
#' object of class \code{METData}.
#'
#' @param METData_training Will be used as training set.
#'
#' @param pheno_new \code{data.frame} object with at least 3 columns.
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'   }
#'   Basically, this `data.frame` only indicates which genotypes should be
#'   predicted in which environments. No columns with phenotypic values is
#'   expected.
#'
#' @param info_environments_to_predict \code{data.frame} object with at least the 4 first
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
#'   package (setting argument `ECs_computed = T` in [create_METData()] and using
#'   subsequently the function [get_ECs()]).}
#'
#'
#' @param scenario_weather_data \code{character}. Options are:
#' * **`manually`**: names of columns in the `env_data_manual` data.frame should be the
#'    same as in METData_training$env_data.
#' * **`use_real_EC`**: environmental variables will be estimated for a set of "past" environments based on
#'  information provided in info_environments_to_predict.
#'    using the [use_real_EC()] function. This assumes that the environmental data can be
#'   retrieved and are from the past. only possible if...
#'   different according to compute_EC.
#' * **`mean_previous_years`**: ECs computed based on average over 5 last years
#' of the location indicated in info_environments_to_predict
#'   use of the location data from the location in the training set (implies
#'   that ) if computed_ECs== TRUE in METData_training.
#'   give fictive harvest date --> example
#' * **`previous_year`**:
#'
#' @param location_to_use \code{character}. Default is NULL.
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
#'   * \strong{Providing env_data  and setting `compute_climatic_ECs = T` is
#'   possible. For instance, the user can have some information regarding the
#'   soil composition (\% % clay, sand, silt, organic matter content).
#'   A disease status can also be encoded as categorical variable if it affected
#'   some environments. In addition to these covariates, weather-based
#'   covariates will be computed if `compute_climatic_ECs = T` with the package if they
#'   were not available to provide as input by the user.}
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
#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#' @examples
#'
#' data(geno_japonica)
#' data(map_japonica)
#' data(pheno_japonica)
#' pheno_japonica2<-pheno_japonica[which(pheno_japonica$year==2013),]
#' data(info_environments_japonica)
#' info_environments_japonica2<-info_environments_japonica[which(info_environments_japonica$year==2013),]
#' data(env_data_japonica)
#' env_data_japonica2<-env_data_japonica[which(env_data_japonica$year==2013),]
#'
#'
add_new_METData <-
  function(METData_training,
           pheno_new = NULL,
           geno_new = NULL,
           compute_climatic_ECs = FALSE,
           info_environments_to_predict = NULL,
           env_data_manual = NULL,
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
    
    
    
    if (compute_climatic_ECs |
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
    
    
    
    if (compute_climatic_ECs) {
      cat('Computation of environmental covariates starts for untested environments.\n')
      merged_ECs <-
        get_ECs(info_environments = info_environments_to_predict,
                ...)
      
      
      # Add ECs to the table env_data, if this table already contains
      # environmental covariates.
      
      if (!is.null(env_data_manual)) {
        env_data_new <-
          merge(merged_ECs,
                env_data_manual %>% select(-location, -year),
                by = 'IDenv')
      }
      
      
      
      else{
        env_data_new <- merged_ECs
      }
      
      ECs_computed <- TRUE
      
      
      cat('Computation of environmental covariates for environments to predict is done.\n')
      
    } else if (!compute_climatic_ECs) {
      ECs_computed <- FALSE
    }
    
    
    ################################
    ## PREPARE GENOMIC DATA ##
    ################################
    
    
    ## Test the geno_new data input ##
    
    if (is.null(geno_new) &
        all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
      geno_new <- METData_training$geno
      map_new <- METData_training$map_markers
      
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
        METData_training$map_markers[which(METData_training$map_markers$marker_name %in%
                                             common_cols), ]
      
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
        info_environments_to_predict %>% dplyr::select(-year, -location, -longitude, -latitude),
        by = 'IDenv',
        all.x = T
      )
      
    }
    
    
    
    
    METData <- structure(
      list(
        'geno' = geno_new,
        'map_markers' = map_new,
        'pheno' = pheno_new,
        'compute_climatic_ECs' = compute_climatic_ECs,
        'ECs_computed' = ECs_computed,
        'env_data' = env_data_new,
        'info_environments' = info_environments_to_predict
      ),
      class = 'METData'
    )
    
    return(METData)
    
    
  }
