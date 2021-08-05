#' @author Cathy C. Jubin \email{cathy.jubin@@uni-goettingen.de}
#' @export
#'
snp_based_haploblocks <- function(geno, map, min_freq = 0.05, k) {
  # convert to matrix
  geno <- as.matrix(geno)
  t_geno <- as.data.frame(t(geno))
  
  map$chr <- as.factor(map$chr)
  map_chr <- split(map, map$chr)
  
  t_geno$chr <- map[match(rownames(t_geno), map$marker_name), 'chr']
  t_geno_chr <- split(t_geno, t_geno$chr)
  t_geno_chr <-
    lapply(t_geno_chr, function(x)
      x[, -which(colnames(x) == 'chr')])
  geno_chr <- lapply(t_geno_chr, t)
  names(geno_chr) <- paste0('chr', unique(map$chr))
  
  
  
  
  haplo_alleles_by_ind_all_chr <- function(x, k) {
    print(x)
    haplotypes_by_chr <- function(nb_chr, geno_chr) {
      geno_data <- geno_chr[[paste0('chr', nb_chr)]]
      genotype_by_window <-
        as.data.frame(zoo::rollapply(geno_data[x, ], k, paste, by = k))
      genotype_by_window$window_nb <-
        paste0('chr', nb_chr, '_W', 1:nrow(genotype_by_window))
      
      
      one2two <- function(x, geno, k) {
        macvec <- as.matrix(geno[x,-(k + 1)])
        #print(class(macvec))
        if (sum(macvec == 1) > 2) {
          a1 <- rep(NA, k)
          a2 <- rep(NA, k)
        } else{
          a1 <-
            ifelse(is.na(macvec[-(k + 1)]), NA, ifelse(macvec[-(k + 1)] == 0, 1, ifelse(macvec[-(k +
                                                                                                   1)] == 2, 2, 1)))
          
          a2 <-
            ifelse(is.na(macvec[-(k + 1)]), NA, ifelse(macvec[-(k + 1)] == 0, 1, ifelse(macvec[-(k +
                                                                                                   1)] == 2, 2, 2)))
        }
        
        haplotypes <- data.table::as.data.table(rbind(a1, a2))
        haplotypes$window_nb <- rep(geno[x,'window_nb'], 2)
        #haplotypes$haplotype <-
        #  paste0(macvec['window_nb'], '_', c('H1', 'H2'))
        haplotypes <-
          haplotypes %>%  unite("haplo_allele", remove = FALSE)
        return(haplotypes)
      }
      
      
      two <-
        lapply(seq_len(nrow(genotype_by_window)), function(x) {
          one2two(x = x,
                  geno = genotype_by_window,
                  k = k)
        })
      twomat <- do.call("rbind", two)
    }
    
    one_ind_all_chr <-
      lapply(seq_len(length(unique(map$chr))), function(x) {
        haplotypes_by_chr(nb_chr = x, geno_chr = geno_chr)
      })
    all_windows_one_ind <- do.call("rbind", one_ind_all_chr)
    
    all_windows_one_ind$ind <- rownames(geno_chr[[1]])[x]
    return(all_windows_one_ind)
  }
  
  all_haplotypes_all_ind <-
    lapply(seq_len(nrow(geno)), function(x) {
      haplo_alleles_by_ind_all_chr(x = x, k = k)
    })
  all_haplotypes <- do.call("rbind", all_haplotypes_all_ind)
  
  freq_table <-
    as.data.frame(table(all_haplotypes$haplo_allele, all_haplotypes$ind))
  
  freq_table <- data.table::as.data.table(freq_table)
  design_mat <-
    data.table::dcast(freq_table, Var2 ~ Var1, value.var = "Freq")
  
  ## Change names of haplotype alleles, and create a table of haplotypes (SNPs,
  ## chromosome and position)
  
  haplo_alleles <-
    colnames(design_mat)[colnames(design_mat) != 'Var2']
  identifier_window <-
    paste0('chr', sub('.*\\_chr', '', haplo_alleles))
  haplo_alleles <- data.frame(cbind(haplo_alleles, identifier_window))
  haplo_alleles$identifier_window <-
    as.factor(haplo_alleles$identifier_window)
  haplo_alleles <-
    as.data.frame(haplo_alleles %>% group_by(identifier_window) %>% mutate(count = row_number(identifier_window)))
  haplo_alleles$newID <-
    paste0(haplo_alleles$identifier_window,
           '_Allele',
           haplo_alleles$count)
  
  colnames(design_mat) <-
    c('geno_ID', haplo_alleles[match(colnames(design_mat)[colnames(design_mat) !=
                                                            'Var2'], haplo_alleles$haplo_alleles), 'newID'])
  
  # Get windows
  names_snps <- list()
  for (s in as.numeric(unique(map$chr))) {
    names_snps[[s]] <-
      as.data.frame(zoo::rollapply(colnames(geno_chr[[s]]), k, paste, by = k))
    names_snps[[s]]$identifier_window <-
      paste0('chr', s, '_W', 1:nrow(names_snps[[s]]))
    
  }
  snps_windows <- do.call("rbind", names_snps)
  haplo_alleles <-
    unique(haplo_alleles[, c('identifier_window', 'newID')])
  haplo_alleles <-
    merge(haplo_alleles, snps_windows, by = 'identifier_window')
  
  ## Remove haplotypes with low variance or with haplotype frequency below a certain level
  
  rownames(design_mat) <- design_mat$Var2
  design_mat <- design_mat[, Var2 := NULL]
  p <- colSums(design_mat) / (2 * nrow(design_mat))
  maf <- 1 - p
  to_keep <- which(maf > min_freq)
  to_keep <- colnames(design_mat)[to_keep]
  design_mat <- design_mat[, ..to_keep]
  
  return(list(
    'haplotype_matrix' = design_mat,
    'info_haplo_alleles' = haplo_alleles
  ))
}