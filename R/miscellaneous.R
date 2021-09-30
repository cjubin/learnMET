#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.
#' @format A \code{data.frame} object with the 4 following columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'  }
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' @name info_environments_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.
#' @format \code{data.frame} object with 7 columns:
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'     \item GY \code{character} Grain Yield
#'     \item PHR \code{character} Percentage of Head Rice Recovery 
#'     \item GC \code{character} Percentage of Chalky Grain
#'     \item PH \code{character} Plant Height
#'   }
#'   
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' @name pheno_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.
#' @format \code{numeric} matrix with genotype values stored in a \code{matrix} or
#'   \code{data.frame} which contains the geno_ID as row.names and markers as
#'   columns.
#'   Dimensions are `n`=327 and `p`=92430.
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' @name geno_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL