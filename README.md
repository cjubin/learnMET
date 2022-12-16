
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **learnMET**

### **Current Version**: 1.1.0 (29 September 2022)

[![DOI](https://img.shields.io/badge/DOI-doi.org%2F10.1093%2Fg3journal%2Fjkac226-B31B1B.svg)](https://doi.org/10.1093/g3journal/jkac226)
[![Release
date](https://img.shields.io/github/release-date/cjubin/learnMET)](https://packagist.org/packages/cjubin/learnMET)
[![GitHub
tag](https://img.shields.io/github/tag/Naereen/StrapDown.js.svg)](https://GitHub.com/Naereen/StrapDown.js/tags/)

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
<https://cjubin.github.io/learnMET/>  
Vignettes are displayed under the Articles section.

# Publication

A publication is available that describes the main features of the
package and how to apply the different functions as a workflow. Results
are provided for several Machine Learning state-of-the-art models tested
with two breeding datasets:  

-   learnMET: an R package to apply machine learning methods for genomic
    prediction using multi-environment trial data Cathy C. Westhues,
    Henner Simianer, Timothy M. Beissinger. G3. doi:
    <https://doi.org/10.1093/g3journal/jkac226>

# Feedback

We are glad about any new user testing learnMET!  
Please contact us if you encounter issues to use some functions of the
package (contact: <cathy.jubin@hotmail.com>).  
Please also do not hesitate to report errors, or additional features
that could be added to the package.
