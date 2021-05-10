% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select_markers.R
\name{select_markers}
\alias{select_markers}
\title{Selection of specific SNPs covariates.}
\usage{
select_markers(
  METData,
  trait,
  method_marker_effects = c("elasticnet", "FarmCPU"),
  method_selection_EN = c("only_variance_across_env", "effect_size_per_env",
    "combined_selection"),
  size_subset_most_variable_markers = 200,
  size_top_markers_by_env = 50,
  plot_penalty_regression_coefficients = F,
  plot_gwas = T,
  path_save_plot = NULL,
  path_save_results = NULL,
  ...
)
}
\arguments{
\item{method_marker_effects}{\code{character} Name of the method to estimate
marker effects in each environment.}

\item{method_selection_EN}{\code{character} Name of the method to select
markers kept for further analyses.}

\item{size_subset_most_variable_markers}{\code{numeric} Number of markers
kept if the selection is based on the variability of marker effects across
environments.}

\item{size_top_markers_by_env}{\code{numeric} Number of markers kept if the
selection is based on marker effect size by environment.}

\item{plot_penalty_regression_coefficients}{\code{logical} Whether to plot
on a grid environment ~ Chromosome the results from the Elastic Net
variable selection by env.}

\item{plot_gwas}{\code{logical} Whether to plot on a grid
environment ~ Chromosome the results from the GWAS by env.}

\item{path_save_plot}{\code{character} Path where the plot should be saved.}

\item{path_save_results}{\code{character} Path where the results (from EN or
from GWAS) should be saved.}

\item{METData.}{An object created by the initial function of the package,
"create_METData.R"}

\item{trait.}{\code{character} Name of the trait under study for which a
subset of markers should be chosen.}
}
\value{
a \code{list} of class \code{METData} which contains the following
elements:
\describe{

\item{geno}{\code{matrix} with genotype values of phenotyped individuals.}

\item{map}{\code{data.frame} with genetic map.}

\item{pheno}{\code{data.frame} with phenotypic trait values.}

\item{compute_EC_by_geno}{\code{Logical} indicates if environmental
covariates should be later computed.}

\item{env_data}{\code{data.frame} with the environmental covariates per
environment (and if genotype-specific, per genotype).}

\item{info_environments}{\code{data.frame} contains basic information on
each environment.}

\item{unique_EC_by_geno}{\code{Logical} to indicate if the EC is genotype-
specific.}

\item{filtering_markers}{\code{Logical} indicates if a filtering marker step
should be applied in further steps}

\item{selected_markers}{\code{character}. Vector containing the names of the
markers selected for further analyses}
}
}
\description{
Select markers based on either:
\enumerate{
\item their effect size and/or on the variance of
their effects across environments estimated by a penalized linear regression
model.
\item GWAS in each environment implemented via FarmCPU
}
}
\author{
Cathy C. Jubin \email{cathy.jubin@uni-goettingen.de}
}