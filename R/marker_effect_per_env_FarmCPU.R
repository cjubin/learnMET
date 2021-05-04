#' Compute marker P-values for each environment with FarmCPU (GWAS method). 
#' Multiple testing correction: Benjamini–Hochberg procedure with alpha=0.05.
#'
#' @param geno \code{data.frame} with markers in columns and inviduals in rows.
#'
#'   Typical input is METData$geno, after having created a METData object using
#'   create_METData() function
#'
#'
#' @param pheno \code{data.frame} with:
#'   First column: ID genotypes
#'   Second column: year
#'   Third column: location
#'   Subsequent columns: phenotypic traits with names indicated in colnames()
#'   Last column: IDenv (combination LocationxYear)
#'
#'   Typical input is METData$pheno, after having created a METData object using
#'   create_METData() function
#'
#' @param map \code{data.frame} with 3 columns.
#'   First column \code{character} with marker names
#'   Second column \code{numeric} with chromosome number
#'   Third column \code{numeric} with marker position.
#'
#'   Typical input is METData$map, after having created a METData object using
#'   create_METData() function
#'
#'
#' @param environment \code{character} indicating the name of the environment
#'   for which marker effects should be computed
#'
#' @param pheno_trait \code{character} Name of the trait under study for which marker
#'   effects should be estimated
#'
#' @param nb_pcs \code{numeric} Number of principal components to use in FarmCPU
#'   analyses.
#'
#' @return a \code{list} which contains the following elements:
#'
#' \describe{
#'   \item{GWAS_results}{\code{data.frame} FarmCPU results for all SNPs}
#'   \item{threshold}{\code{numeric} Cutoff derived from Benjamini-Hochberg
#'   procedure}
#'   \item{selected_markers}{\code{character} Vector containing the names of the
#'   SNPs passing the threshold}
#' }


marker_effect_per_env_FarmCPU <-
  
  function(geno,
           pheno,
           map,
           environment,
           pheno_trait,
           nb_pcs = 5,
           ...) {
    # Select the phenotype data corresponding to the selected environment
    
    pheno <- pheno[pheno$IDenv == environment, ]
    
    list_predictors <- colnames(geno)
    
    geno$geno_ID = row.names(geno)
    pheno <- merge(pheno, geno, by = 'geno_ID', all.x = T)
    
    # Select trait and marker columns from the phenotypic file merged with geno
    # data.
    
    pheno <- pheno[, c(pheno_trait, list_predictors, 'geno_ID')]
    colnames(pheno)[ncol(pheno)] <- 'Taxa'
    colnames(geno)[ncol(geno)] <- 'Taxa'
    
    
    colnames(map) <- c('SNP', 'Chromosome', 'Position')
    map <- map[map$SNP %in% colnames(geno),]
    
    # Association analysis in an environment
    
    myGAPIT <- GAPIT3::GAPIT(
      Y = pheno[, c('Taxa', pheno_trait)],
      GD = pheno[, c('Taxa', colnames(geno)[colnames(geno) %notin% 'Taxa'])],
      GM = map,
      model = c("FarmCPU"),
      Geno.View.output = FALSE,
      PCA.total = nb_pcs,
      file.output = F
    )
    
    # Benjamini Hochberg procedure
    # Order the p-values
    
    pValues <- sort(myGAPIT$GWAS$P.value)
    
    # BH-procedure with alpha = 0.05
    alpha <- 0.05
    m <- length(pValues)
    qValues <- c()
    
    for (i in 1:m) {
      qV <- (i / m) * alpha
      qValues <- append(qValues, qV)
    }
    
    # find the largest p-value satisfying p_i < q_i
    
    BH_test <- qValues > pValues
    
    threshold <- pValues[sum(BH_test)]
    
    selected_markers <-
      myGAPIT$GWAS[which(myGAPIT$GWAS$P.value <= threshold), 'SNP']
    
    gwas_results <- myGAPIT$GWAS
    
    gwas_results$environment <- environment
    
    if (length(threshold) == 0) {
      gwas_results$threshold <- NA
    } else{
      gwas_results$threshold <- threshold
    }
    
    return(
      list(
        'GWAS_results' = gwas_results,
        'threshold' = threshold,
        'selected_markers' = selected_markers
      )
    )
    
  }
