#' Select marker effects based on either:
#' 1. their effect size and/or on the variance of
#' their effects across environments estimated by a penalized linear regression
#' model.
#' 2. GWAS in each environment implemented via FarmCPU
#'
#' @param METData. An object created by the initial function of the package,
#'   "create_METData.R"
#'
#' @param trait. \code{character} Name of the trait under study for which a
#'   subset of markers should be chosen.
#'
#' @param method_marker_effects \code{character} Name of the method to estimate
#'   marker effects in each environment.
#'
#' @param method_selection_EN \code{character} Name of the method to select 
#'   markers kept for further analyses.
#'
#' @param size_subset_most_variable_markers \code{numeric} Number of markers
#'  kept if the selection is based on the variability of marker effects across 
#'  environments.
#'
#' @param size_top_markers_by_env \code{numeric} Number of markers kept if the
#'   selection is based on marker effect size by environment.
#' 
#' @param plot_penalty_regression_coefficients \code{logical} Whether to plot
#'   on a grid environment ~ Chromosome the results from the Elastic Net 
#'   variable selection by env.
#' 
#' @param plot_gwas \code{logical} Whether to plot on a grid 
#'   environment ~ Chromosome the results from the GWAS by env.
#'  
#' @param path_save_plot \code{character} Path where the plot should be saved.
#'   
#' @param path_save_results \code{character} Path where the results (from EN or
#' from GWAS) should be saved.
#'
#' @return a \code{list} of class \code{METData} which contains the following 
#' elements:
#' \describe{
#' 
#'   \item{geno}{\code{matrix} with genotype values of phenotyped individuals.}
#'
#'   \item{map}{\code{data.frame} with genetic map.}
#'
#'   \item{pheno}{\code{data.frame} with phenotypic trait values.}
#'
#'   \item{compute_EC_by_geno}{\code{Logical} indicates if environmental
#'   covariates should be later computed.}
#'
#'   \item{env_data}{\code{data.frame} with the environmental covariates per
#'   environment (and if genotype-specific, per genotype).}
#'
#'   \item{info_environments}{\code{data.frame} contains basic information on
#'   each environment.}
#'
#'   \item{unique_EC_by_geno}{\code{Logical} to indicate if the EC is genotype-
#'   specific.}
#'
#'   \item{filtering_markers}{\code{Logical} indicates if a filtering marker step
#'   should be applied in further steps}
#'
#'   \item{selected_markers}{\code{character}. Vector containing the names of the
#'   markers selected for further analyses}
#' }



select_markers <- function(METData,
                           trait,
                           method_marker_effects = c('elasticnet', 'FarmCPU'),
                           method_selection_EN = c('only_variance_across_env',
                                                'effect_size_per_env',
                                                'combined_selection'),
                           size_subset_most_variable_markers = 200,
                           size_top_markers_by_env = 50,
                           plot_penalty_regression_coefficients = F,
                           plot_gwas = T,
                           path_save_plot = NULL,
                           path_save_results = NULL,
                           ...) {
  # If the number of markers is less than 1000, all markers can be used
  # in subsequent analyses
 
    
  if (dim(METData$geno)[2] < 1000) {
    stop('The number of markers is low and does not need to be further reduced.')
  }
  
  
  # Vector containing names of all environments in the MET analysis
  
  all_envs = unique(METData$info_environments$IDenv)
  
  
  # According to the selected method for calculating marker effects with CV,
  # the marker effects in each environment are computed.
  
  if (method_marker_effects == 'elasticnet') {
    res_all_envs = lapply(
      all_envs,
      function (x) {marker_effect_per_env_EN(environment = x,
                                             geno = METData$geno,
                                             pheno = METData$pheno,
                                             pheno_trait = trait,
                                             ...)})
    
    
    saveRDS(res_all_envs,file=paste0(path_save_results,'/elasticnet_results.RDS'))
    
    # Method 1 for EN: evaluate the variance of marker effects across environments and
    # select accordingly a subset of markers of a certain size.
    
    if (method_selection_EN == 'only_variance_across_env') {
      marker_effects_all_env <- do.call("rbind", res_all_envs)
      
      variance_markers_across_env <-
        marker_effects_all_env %>% group_by(term) %>%
        summarise(var = var(cv_mean))
      
      selected_markers <-
        top_n(variance_markers_across_env,
              size_subset_most_variable_markers,
              var)[, 1]
      
      selected_markers <- selected_markers$term
      
    }
    
    # Method 2 for EN: select the top n markers the most significant in each environment.
    
    if (method_selection_EN == 'effect_size_per_env') {
      marker_effects_all_env <- do.call("rbind", res_all_envs)
      
      subset_top_markers_by_env <-
        marker_effects_all_env %>% group_by(environment) %>%
        top_n(size_top_markers_by_env,
              cv_mean)
      
      selected_markers <- unique(subset_top_markers_by_env[, 1])
      
      selected_markers <- selected_markers$term
      
    }
    
    if (plot_penalty_regression_coefficients == T) {
      # Merge marker name with the table of marker positions + chromosome
      
      markers_table <-
        merge(
          marker_effects_all_env,
          METData$map_markers,
          by.x = 'term',
          by.y = 'marker_name'
        )
      
      # Grid plot with chromosome along columns and environments in rows
      
      plot1 <-
        ggplot(data = markers_table, aes(x = pos, y = cv_mean)) +
        geom_point() +
        theme_bw() +
        facet_grid(environment ~ chr) +
        geom_text_repel(
          data = subset(markers_table, cv_mean > 0.01),
          aes(label = term),
          size = 4
        ) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        labs(title = "Marker effects estimated with Elastic Net within each environment by CV", x =
               "SNP position", y = "Penalty regression coefficient")
      
      print(plot1)
      
      ggsave(plot1,
             filename = paste0(path_save_plot, '/plot_marker_effecs_EN_by_env.pdf'))
      
    } 
  }
  
  if (method_marker_effects == 'FarmCPU') {
    
    res_all_envs = lapply(
      all_envs,
      function (x) {marker_effect_per_env_FarmCPU(environment = x,
                                               geno = METData$geno,
                                               pheno = METData$pheno,
                                               pheno_trait = trait,
                                               ...)})
    
    saveRDS(res_all_envs,file=paste0(path_save_results,'/FarmCPU_results.RDS'))
    
    # Select the third element of each list element from res_all_envs:
    # Contains the markers passing the B-H procedure cutoff
    
    selected_markers <- unlist(res_all_envs %>%  map(c('selected_markers')))
    
    # Plot of GWAS results per environment:
    
    if (plot_gwas == T){
      
      # Select the gwas tables from each sub-element of the res_all_envs list
      # Bind the data.frames
      
      gwas_tables <- res_all_envs %>%  map(c('GWAS_results'))
      
      markers_table <- do.call("rbind", gwas_tables)
      
      # Grid plot with chromosome along columns and environments in rows
      
      plot1 <-
        ggplot(data = markers_table, aes(x = Position, y = -log10(P.value))) +
        geom_point() +
        theme_bw() +
        facet_grid(environment ~ Chromosome) +
        geom_text_repel(
          data = subset(markers_table, -log10(P.value)>=-log10(threshold) ),
          aes(label = SNP),
          size = 4
        ) +
        geom_hline(aes(yintercept=-log10(threshold)),linetype='dashed', col = 'red')+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        labs(title = "GWAS results with FarmCPU by environment", x =
               "SNP position", y = "-log10(P.value)")
      
      print(plot1)
      
      ggsave(plot1,
             filename = paste0(path_save_plot, '/plot_gwas_by_env.pdf'))
    }
    
  }
  
  METData <- list(
    'geno' = METData$geno,
    'map_markers' = METData$map,
    'pheno' = METData$pheno,
    'compute_ECs' = METData$compute_ECs,
    'env_data' = METData$env_data,
    'info_environments' = METData$info_environments,
    'unique_EC_by_geno' = METData$unique_EC_by_geno,
    'filtering_markers' = METData$filtering_markers,
    'selected_markers' = selected_markers
  )
  
  class(METData) <- c("METData", "list")
  
  return(METData)
  
}
