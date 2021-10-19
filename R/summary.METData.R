#' Summary of an object of class METData
#' 
#' @description Different basic information about the MET
#' dataset are retrieved from the METData object.
#' @param object of class `METData`
#' @return object of class `summary.METData`
#' @method summary METData
#' @export
summary.METData <- function(object, ...) {
  stopifnot(inherits(object, "METData"))
  
  ans <- list()
  
  # summary of info environments
  
  ans$env$n <- length(unique(object$info_environments$IDenv))
  ans$env$n_unique_years <-
    length(unique(object$info_environments$year))
  ans$env$n_unique_locations <-
    length(unique(object$info_environments$location))
  
  # summary of allocation of phenotypic obs.
  
  ans$pheno$distribution_pheno <-
    table(object$pheno$year, object$pheno$location)
  ans$pheno$unique_geno <- length(unique(object$pheno$geno_ID))
  
  
  # summary of climate info
  ans$climate_data$climate_data_retrieved <-
    object$climate_data_retrieved
  ans$climate_data$list_climatic_predictors <-
    length(object$list_climatic_predictors)
  
  # summary of soil info
  ans$soil_data$list_climatic_predictors <-
    length(object$list_soil_predictors)
  
  # summary of pheno info
  ans$pheno$n_traits <- dim(object$pheno)[[2]] - 4
  ans$pheno$stats <-
    summary(data.frame(apply(
      object$pheno %>% dplyr::select(-geno_ID, -year, -location,-IDenv),
      2,
      as.numeric
    )))
  
  # summary geno matrix
  ans$geno$p <- ncol(object$geno)
  
  # summary map
  
  if (!is.null(object$map)) {
    ans$map$nb_chr <- length(table(object$map$chr))
    ans$map$table <- table(x$map$chr)
  }
  
  class(ans) <- "summary.METData"
  return(ans)
}


#' Print the summary of an object of class METData
#' 
#' @description Print different basic information about the MET
#' dataset are retrieved from the METData object.
#' @param object of class `summary.METData``
#' @method print summary.METData
#' @export

print.summary.METData <- function(x, ...) {
  cat("x of class 'METData' \n")
  cat("--------------------------\n")
  cat("General information about the MET data \n")
  cat("\n")
  cat("No. of unique environments represented in the data:",
      x$env$n,
      "\n")
  cat("       unique years represented in the data:",
      x$env$n_unique_years,
      "\n")
  cat("       unique locations represented in the data:",
      x$env$n_unique_locations,
      "\n")
  cat("\n")
  
  cat("Distribution of phenotypic observations according to year and location: \n")
  print(x$pheno$distribution_pheno)
  cat("No. of unique genotypes which are phenotyped\n",
      x$pheno$unique_geno,
      "\n")
  
  cat("--------------------------\n")
  cat("Climate variables \n")
  cat("Weather data extracted from NASAPOWER?\n")
  if (x$climate_data$climate_data_retrieved) {
    cat('YES\n')
  } else{
    cat('NO\n')
  }
  cat(
    "No. of climate variables available:",
    x$climate_data$list_climatic_predictors,
    "\n"
  )
  
  
  cat("--------------------------\n")
  cat("Soil variables \n")
  
  cat("No. of soil variables available:",
      x$soil_data$list_climatic_predictors,
      "\n")
  
  
  
  cat("--------------------------\n")
  
  cat("Phenotypic data \n")
  cat("\t No. of traits: ", x$pheno$n_traits, "\n")
  cat("\n")
  print(x$pheno$stats)
  
  
  
  cat("--------------------------\n")
  
  cat("Genotypic data \n")
  cat("\t No. of markers :", x$geno$p, "\n")
  
  if (!is.null(x$map)) {
    cat("--------------------------\n")
    cat("Map data\n")
    cat("\t No. of chromosomes    ", x$map$nb_chr, "\n\n")
    cat("\t markers per chromosome \n\t")
    print(x$map$table)
  }
  invisible(x)
  
  
}
