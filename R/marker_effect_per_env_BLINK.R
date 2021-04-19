
#' Compute marker effects per environment with BLINK
#'
#'
#' @param environment \code{character} indicating the name of the environment
#' for which marker effects should be computed
#'
#' @param geno \code{data.frame} with markers in columns and inviduals in rows.
#' Typical input is METData$geno.
#'
#' @param pheno \code{data.frame} with:
#' First column: ID genotypes
#' Second column: year
#' Third column: location
#' Subsequent columns: phenotypic traits with names indicated in colnames()
#' Last column: IDenv (combination LocationxYear)
#' Typical input is METData$pheno
#'
#' @param trait \code{character} Name of the trait under study for which marker
#' effects should be estimated
#'
#'
#' @return 
#' First column contains the marker names.
#' Second column contains the marker effects in this environment
#' calculated by cross-validation.
#' Third column contains the environment name (combination LocationxYear).
#'


marker_effect_per_env_BLINK <-
  
  function(environment,
           geno,
           pheno,
           map,
           trait,
           nb_folds_cv = 10,
           reps = 2) {
    
    # Select the phenotype data corresponding to the selected environment
    
    pheno <- pheno[pheno$IDenv == environment, ]
    list_predictors <- colnames(geno)
    geno$geno_ID = row.names(geno)
    pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)
    
    # Select trait and marker columns from the phenotypic file merged with geno
    # data.
    
    pheno <- pheno[, c(trait, list_predictors,'geno_ID')]
    colnames(pheno)[ncol(pheno)] <- 'Taxa'
    colnames(geno)[ncol(geno)] <- 'Taxa'
    
    # Association analysis in an environment
    
    colnames(map) <- c('SNP','Chromosome','Position')
    map$Chromosome <- as.numeric(map$Chromosome)
    
    myGAPIT_MLM <- GAPIT3::GAPIT(
      Y=pheno[,c('Taxa',trait)],
      GD=pheno[,c('Taxa',colnames(geno)[colnames(geno)%notin%'Taxa'])],
      GM=map,
      model=c("SUPER"),
      Geno.View.output=FALSE,
      #PCA.total=3
    )
    
    return(gwas_results)
    
  }


