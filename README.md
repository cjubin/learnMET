
<!-- README.md is generated from README.Rmd. Please edit that file -->

# learnMET

[![Github All
Releases](https://img.shields.io/github/downloads/cjubin/learnMET/total.svg)]()

`learnMET` (**learn** **M**ulti-**E**nvironment **T**rials) provides a
pipeline for crop predictive breeding. In particular, `learnMET` (1)
facilitate environmental characterization via the retrieval and
aggregation of daily weather data; (2) allows the evaluation of various
types of state-of-the-art machine learning approaches based on relevant
cross-validation schemes for multi-environment trial datasets (3)
enables to implement predictions for unobserved configurations of
genotypic and environmental predictors that the user wants to test *in
silico*.  
In the Reference section, the different functions implemented in the
package are listed. **Only the so called main functions have to be run
by the user in a typical workflow**.

# Installation

Install the development version from
[GitHub](https://github.com/cjubin/learnMET) with:

``` r
devtools::install_github("cjubin/learnMET")

# To build the HTML vignette use
devtools::install_github("cjubin/learnMET", build_vignettes = TRUE)
```

# Package documentation and vignettes

Vignettes and documentation are available at:
<https://cjubin.github.io/learnMET/> Vignettes are displayed under the
Articles section.

# Preprint

A preprint is available that describes the main current features of the
package: \* learnMET: an R package to apply machine learning methods for
genomic prediction using multi-environment trial data Cathy C. Westhues,
Henner Simianer, Timothy M. Beissinger. bioRxiv 2021.12.13.472185; doi:
<https://doi.org/10.1101/2021.12.13.472185>

# Feedback

We are glad about any new user testing learnMET! Please contact us if
you encounter issues to create a METData object based on your data, or
to run a prediction model, etc. Please also do not hesitate to report
errors, or additional features that could be added to the package.
