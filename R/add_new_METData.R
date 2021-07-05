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
#'   package (setting argument `ECs_computed = T` in [create_METData()] and using
#'   subsequently the function [get_ECs()]).}
#'
#'
#' @param unique_EC_by_geno \code{logical} Default is FALSE. TRUE not yet
#'   implemented.
#'
#' @param scenario_weather_data \code{character}. Options are:
#' * **`manually`**: names of columns in the `env_data_manual` data.frame should be the
#'    same as in original METData$env_data).
#' * **`use_real_EC`**: environmental variables will be estimated for a set of "past" environments based on
#'  information provided in info_environments_new.
#'    using the [use_real_EC()] function. This assumes that the environmental data can be
#'   retrieved and are from the past. only possible if...
#'   different according to compute_EC.
#' * **`mean_previous_years`**: ECs computed based on average over 5 last years
#' of the location indicated in info_environments_new
#'   use of the location data from the location in the training set (implies
#'   that ) if computed_ECs== TRUE in METData_training.
#'   give fictive harvest date --> example
#' * **`previous_year`**: 
#'
#' @param location_to_use \code{character}. Default is NULL.
#'
#' @param env_data_manual \code{data.frame} 
#' can contain soil variables (with same column names as used in METData_training$env_data) but no weather-based covariates
#' 
#' can be let as NULL by user, if no
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

#'
#' @param unique_EC_by_geno \code{logical} indicates if the environmental
#'   covariates contained in env_data_manual are also genotype-specific (dependent on
#'   the phenology, for instance) or unique for a complete environment. Default
#'   is `FALSE`.
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
           geno_new = NULL,
           pheno_new,
           info_environments_new = NULL,
           unique_EC_by_geno = FALSE,
           scenario_weather_data,
           env_data_manual = NULL,
           crop = NULL) {
    
    ## Test that the pheno data to predict are provided ##
    
    if (is.null(pheno_new)) {
      stop("Pheno data to predict not provided via pheno_new argument")
    }
    
    if (scenario_weather_data == 'use_real_EC' & is.null(info_environments_new)) {
      stop(paste("Information about environments in info_environments_new",
                 "needed to derive environmental covariates with argument",
                 "scenario_weather_data ==  'use_real_EC"))
    } else if (scenario_weather_data == 'use_real_EC' & !is.null(info_environments_new) & ncol(info_environments_new)<6){
      stop('All columns must be indicated in info_environments_new')
    }
    
    if (scenario_weather_data %in% c('mean_previous_years','previous_year') & METData_training$compute_ECs & is.null(info_environments_new)) {
      stop(paste("Information about environments in info_environments_new",
                 "needed to derive environmental covariates with argument",
                 "scenario_weather_data ==  'use_real_EC"))
    } else if (scenario_weather_data %in% c('mean_previous_years','previous_year') & METData_training$ECs_computed & !is.null(info_environments_new) & ncol(info_environments_new)<6){
      stop('Longitude and latitude data for the locations must be provided in info_environments_new.')
    }
    
    if (scenario_weather_data != 'use_real_EC' & is.null(info_environments_new)){
      info_environments <- unique(pheno[,c('location','year')])
      info_environments$IDenv <-
        paste0(info_environments$location, '_', info_environments$year)
      
    }

    
    ## Test the geno_new data input ##
    
    if (is.null(geno_new) &
        all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
      geno <- METData_training$geno
    } else if (is.null(geno_new) &
               !all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
      stop(
        'Some genotype names are not present in the row.names of METData_training$geno, and no additional genotype data were provided.
           Please add geno data for the lines to be predicted which are not in row.names(METData_training$geno)'
      )} else if (
      !is.null(geno_new) &
        !all(pheno_new$geno_ID %in% row.names(METData_training$geno))) {
          common_cols <- intersect(colnames(METData_training$geno), colnames(geno_new))
          geno <- unique(rbind(METData_training$geno[, common_cols],
                            geno_new[, common_cols]))
          map <- METData_training$map_markers[METData_training$map_markers$marker_name%in%common_cols,] } else {
      geno <- METData_training$geno
      map <- METData_training$map_markers
      
    }
    
    ## Test the pheno new data imput ##
    
    # test correct class for the different columns of the phenotype data
    if (!is.data.frame(pheno_new)) {
      
    }
    if (ncol(pheno_new) < 3) {
      stop(
        'pheno_new_data should contain 3 columns: genotype lines (col1), year (col2), location (col3).'
      )
    }
    
    if (!is.character(pheno_new[, 1])) {
      stop("the genotype names/IDs (first column of pheno_new) must be character")
    }
    if (!is.numeric(pheno_new[, 2])) {
      stop("the year (second column of pheno_new) must be numeric")
    }
    if (!is.character(pheno_new[, 3])) {
      stop("the location (third column of pheno_new) must be character")
    }
    
    
    # Assign col.names to pheno_new columns and transform year + location to 
    # factor
    
    colnames(pheno_new)[1:3] <- c('geno_ID', 'year', 'location')
    pheno_new$year = as.factor(pheno_new$year)
    pheno_new$location = as.factor(pheno_new$location)
    
    
    # Bind with the training METData
    pheno <- plyr::rbind.fill(METData_training$pheno,pheno_new)
    
    # Create unique ID environment based on the location x year combination
    pheno$IDenv <- paste0(pheno$location, '_', pheno$year)
    
    
    
    ## Determine the environmental data based on the chosen scenario_weather_data
    ## option.
    
    if (scenario_weather_data == 'get_EC') {
      cat('Computation of environmental covariates starts.\n')
      env_data <- get_ECs(unique_EC_by_geno = unique_EC_by_geno,
                          env_data = env_data_manual,
                          info_environments = info_environments_new,
                          crop = crop,
                          ...)
      
      ECs_computed <- TRUE
      cat('Computation of environmental covariates is done.\n')
      common_cols <- intersect(colnames(METData_training$env_data), colnames(env_data))
      env_data <- unique(rbind(METData_training$env_data[, common_cols],
                               env_data[, common_cols]))
    }
    
    
    
    if (scenario_weather_data == 'manual' & is.null(env_data_manual)){
      stop(paste("Environmental data should be provided by the user with",
                 "argument scenario_weather_data ==  'manual"))
      
    } else if (scenario_weather_data == 'manual'){
      common_cols <- intersect(colnames(METData_training$env_data), colnames(env_data_manual))
      env_data <- unique(rbind(METData_training$env_data[, common_cols],
                           env_data_manual[, common_cols]))
      
    }
    
    
    if (scenario_weather_data == 'mean_previous_years' & METData_training$ECs_computed){
      table_env_previous_years <- info_environments_new[]
      mean_ECs_previous_years(,...)
      
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
    
      
    # if marker matrix is given, test that the marker names are the same in the map and in the marker genotype matrices
    
    if (!is.null(map) & !identical(colnames(geno), map[, 1])) {
      stop("marker names in genotypic data and in map are not the same")
      
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
    if (ECs_computed &
        is.null(info_environments$harvest.date)) {
      stop('Computation of ECs is required but no date for the harvest date.')
    }
    if (ECs_computed &
        is.null(info_environments$planting.date)) {
      stop('Computation of ECs is required but no date for the planting date.')
    }
    
    if (ECs_computed &
        !inherits(info_environments_G2F$harvest.date, 'Date')) {
      stop('planting date in info_environments as Date (format y-m-d).')
    }
    if (ECs_computed &
        !inherits(info_environments_G2F$planting.date, 'Date')) {
      stop('harvest date in info_environments as Date (format y-m-d).')
    }
    
    
    
    
    
    
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
      
      
      env_data <- merge(
        env_data,
        info_environments %>% dplyr::select(-year,-location,-longitude,-latitude),
        by = 'IDenv',
        all.x = T
      )
      
    } else{
      cat('No environmental covariates provided.\n')
    }
    
    if (ECs_computed) {
      cat('Environmental covariates should be determined.\n')
    }
    
    if (!ECs_computed & is.null(env_data)) {
      cat(
        paste(
          'No environmental covariates will be computed nor used using',
          'the package. To allow calculation of ECs, please use the',
          'argument ECs_computed = T.\n'
        )
      )
    }
    
    
    
    METData <- structure(
      list(
        'geno' = geno,
        'map_markers' = map,
        'pheno' = pheno,
        'ECs_computed' = ECs_computed,
        'env_data' = env_data,
        'info_environments' = info_environments,
        'unique_EC_by_geno' = unique_EC_by_geno,
        'filtering_markers' = filtering_markers
      ),
      class = 'METData'
    )
    
    return(METData)
    
    
  }
