#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Information about environments where *japonica* rice population was grown.
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
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name info_environments_japonica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Climate variables data for the *japonica* population grown in one location (Treinta y Tres, Uruguay) for three years.
#' @format \code{data.frame} containing as many rows as info_environments
#'   Columns are:
#'   \enumerate{
#'     \item geno_ID \code{character} with genotype ID
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'     \item Col 3 to Col 56: climatic covariates derived for each crop
#'     growing stage. More explanation on the ECs in the original publication
#'     from Monteverde et al. (2019)
#'   }
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name climate_variables_japonica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Phenotypic data for the *japonica* rice population for four traits.
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
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name pheno_japonica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.
#' @format \code{numeric} matrix with genotype values stored in a \code{matrix} or
#'   \code{data.frame} for the *japonica* rice population. Names of lines (geno_ID) are given as row.names and markers as
#'   columns.
#'   Dimensions are `n`=321 and `p`=44,598.
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name geno_japonica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Physical map (GBS dataset) for the *japonica* population (44,598 markers)
#' @format  map \code{data.frame} object with 3 columns.
#'   \enumerate{
#'   \item marker_name \code{character} with marker names
#'   \item chr \code{numeric} with chromosome number
#'   \item pos \code{numeric} with marker position.
#'   }
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name map_japonica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL


#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Information about environments where *indica* rice population was grown.
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
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name info_environments_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Climate variables data for the *indica* population grown in one location (Treinta y Tres, Uruguay) for three years.
#' @format \code{data.frame} containing as many rows as info_environments
#'   Columns are:
#'   \enumerate{
#'     \item geno_ID \code{character} with genotype ID
#'     \item year \code{numeric} with the year label
#'     \item location \code{character} with the location character
#'     \item Col 3 to Col 56: climatic covariates derived for each crop
#'     growing stage. More explanation on the ECs in the original publication
#'     from Monteverde et al. (2019)
#'   }
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name climate_variables_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Phenotypic data for the *indica* rice population for four traits.
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
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name pheno_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.
#' @format \code{numeric} matrix with genotype values stored in a \code{matrix} or
#'   \code{data.frame} for the *indica* rice population. Names of lines (geno_ID) are given as row.names and markers as
#'   columns.
#'   Dimensions are `n`=327 and `p`=92,430.
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name geno_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Multi-year trial data of rice
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package.\cr
#' Physical map (GBS dataset) for the *indica* population (92,430 markers)
#' @format  map \code{data.frame} object with 3 columns.
#'   \enumerate{
#'   \item marker_name \code{character} with marker names
#'   \item chr \code{numeric} with chromosome number
#'   \item pos \code{numeric} with marker position.
#'   }
#' @md
#' @source Data from the INIA's Rice Breeding Program (Uruguay). Data downloaded from publication Monteverde et al. (2019).
#' @references 
#' \insertRef{monteverde2019integrating}{learnMET} 
#' Instituto Nacional de Investigación Agropecuaria (INIA). Programa Arroz. Estación Experimental INIA Treinta y Tres, Uruguay
#' @name map_indica
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL

#' @title Maize experimental multi-environment data sets (Genomes to Fields Initiative)
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package. \cr
#' Field information about 22 environments from the Genomes to Fields project.
#' @format A \code{data.frame} object with the 4 following columns. \cr
#'   \enumerate{
#'     \item year: \code{numeric} Year label of the environment
#'     \item location: \code{character} Name of the location
#'     \item longitude: \code{numeric} longitude of the environment
#'     \item latitude: \code{numeric} latitude of the environment
#'     \item planting.date: \code{Date} Planting date, format "YYYY-MM-DD"
#'     \item harvest.date: \code{Date} Harvest date, format "YYYY-MM-DD"
#'  }
#' @md
#' @source Subset of the publicly available datasets from the Genomes to Fields Initiative. \cr 
#' All data are publicly available at: https://www.genomes2fields.org/resources/,
#' @references 
#' \insertRef{mcfarland2020maize}{learnMET} 
#' \insertRef{alkhalifah2018maize}{learnMET} 
#' @name info_environments_G2F
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL



#' @title Maize experimental multi-environment data sets (Genomes to Fields Initiative)
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package. \cr
#' Phenotypic data for 22 environments from the Genomes to Fields project.
#' @format \code{data.frame} object with 6 columns:
#'   \enumerate{
#'     \item geno_ID \code{character} contains the genotype identifiers.
#'     \item year \code{numeric} contains the year of the observation.
#'     \item location \code{character} contains the name of the location.
#'     \item pltht \code{character} plant height (cm)
#'     \item yld_bu_ac \code{character} grain yield (bushels per acre)
#'     \item earht \code{character} ear height (cm)
#'   }
#' @md
#' @source Subset of the publicly available datasets from the Genomes to Fields Initiative. \cr 
#' All data are publicly available at: https://www.genomes2fields.org/resources/,
#' @references 
#' \insertRef{mcfarland2020maize}{learnMET} 
#' \insertRef{alkhalifah2018maize}{learnMET} 
#' @name pheno_G2F
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL




#' @title Maize experimental multi-environment data sets (Genomes to Fields Initiative)
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package. \cr
#' Physical map (SNP information generated using a genotyping-by-sequence (GBS) method for the inbred lines) for
#' 106,414 markers.
#' @format \code{data.frame} object with 3 columns.
#'   \enumerate{
#'   \item marker_name \code{character} with marker names
#'   \item chr \code{numeric} with chromosome number
#'   \item pos \code{numeric} with marker position.
#'   }
#' @md
#' @source Subset of the publicly available datasets from the Genomes to Fields Initiative. \cr 
#' All data are publicly available at: https://www.genomes2fields.org/resources/,
#' @references 
#' \insertRef{mcfarland2020maize}{learnMET} 
#' \insertRef{alkhalifah2018maize}{learnMET} 
#' @name map_G2F
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL



#' @title Maize experimental multi-environment data sets (Genomes to Fields Initiative)
#' @description This data are used as toy data for several functions of
#' \pkg{learnMET} package. \cr
#' Physical map (SNP information generated using a genotyping-by-sequence (GBS) method for the inbred lines) for
#' 106,414 markers.
#' @format \code{data.frame} with genotype values of maize hybrids (*in silico* genotype matrix based on inbred lines) 
#' stored, with geno_ID as row.names and markers in columns. 
#' Dimensions are `n`=1,913 and `p`=106,414.
#' @md
#' @source Subset of the publicly available datasets from the Genomes to Fields Initiative. \cr 
#' All data are publicly available at: https://www.genomes2fields.org/resources/,
#' @references 
#' \insertRef{mcfarland2020maize}{learnMET} 
#' \insertRef{alkhalifah2018maize}{learnMET} 
#' @name geno_G2F
#' @docType data
#' @author Cathy C. Westhues \email{cathy.jubin@uni-goettingen.de}
#' @keywords data
NULL