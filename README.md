# Bayesian inference for Markov Branching Process estimation

This R package makes the estimation of Continuous Time Markov Branching processes easy. By taking advantage of the approximate normality of branching processes with large starting populations, accurate parameter inference is possible for arbitrarily complex multi-type processes. In addition, the package allows for branching process parameters to have arbitrary functional relationships with any number of environmental variables. This makes it easy to infer the relationship between cell kinetics and drug doses in pharmacodynamics experiments, for example.

## Getting Started

To install `bpinference`, first install the R `devtools` package:
```
install.packages(devtools)
```
Then install the `bpinference` pacakge:
```
devtools::install_github('jproney/bpinference')
```
This will install `bpinference` and all of its dependencies.

## Package Structure and Documentation
The core source code of `bpinference` is located in the `/R` subdirectory of the repository.
`/man` contains documentation for individual functions in the package. Look here to understand the return values, parameters, and intended use of specific functions in the package.
`/vignettes` contains long-form, annoated documentation for the package. Each vignette contains a `.Rmd` file which will walk you through the steps of setting up, simulating, and estimating a particular branching process model. Look here for concrete examples with detailed explanations.
`/examples` contains more example uses of the package, but with less thorough explanations than the examples in `/vignettes`

## Getting Started
To become familiar with the workflow of the package, have a look at `/vignettes/first-model.Rmd`. Happy estimation!
