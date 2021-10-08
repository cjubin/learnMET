
<!-- README.md is generated from README.Rmd. Please edit that file -->

# learnMET

<!-- badges: start -->
<!-- badges: end -->

`learnMET` (**learn** **M**ulti-**E**nvironment **T**rials) provides an
integrated pipeline for crop predictive breeding. In particular,
`learnMET` (1) facilitate environmental characterization via the
retrieval and aggregation of daily weather data; (2) allows the
evaluation of various types of state-of-the-art machine learning
approaches based on relevant cross-validation schemes for
multi-environment trial datasets (3) enables to implement predictions
for unobserved configurations of genotypic and environmental predictors
that the user wants to test *in silico*.  
In the Reference section, the different functions implemented in the
package are listed. **Only the so called main functions have to be run
by the user in a typical workflow**. Other listed functions are
functions which are called by the main function according to the
parameters specified by the user.

# Installation

Install the development version from
[GitHub](https://github.com/cjubin/learnMET) with:

``` r
devtools::install_github("cjubin/learnMET")

# To build the HTML vignette use
devtools::install_github("cjubin/learnMET", build_vignettes = TRUE)
```
