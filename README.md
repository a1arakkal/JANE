
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JANE

<!-- badges: start -->

[![R-CMD-check](https://github.com/a1arakkal/JANE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a1arakkal/JANE/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/a1arakkal/JANE/branch/master/graph/badge.svg)](https://app.codecov.io/gh/a1arakkal/JANE?branch=master)

<!-- badges: end -->

JANE is a package to fit and simulate latent space network cluster
models using an expectation-maximization algorithm. Enables flexible
modeling of unweighted network data, supporting both directed and
undirected networks, with or without degree heterogeneity. Designed to
handle large networks efficiently, it allows users to explore network
structure through latent space representations, identify clusters within
network data, and simulate models with varying clustering and
connectivity patterns.

## Installation

``` r
# Anticipated release from CRAN
# install.packages("JANE")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("a1arakkal/JANE")
```

## Basic usage

### Simulate a network

``` r
library(JANE)
mus <- matrix(c(-1,-1,1,-1,1,1), 
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)
omegas <- array(c(diag(rep(7,2)),
                  diag(rep(7,2)), 
                  diag(rep(7,2))), 
                  dim = c(2,2,3))
p <- rep(1/3, 3)
beta0 <- 1.0
sim_data <- JANE::sim_A(N = 100L, 
                        model = "NDH",
                        mus = mus, 
                        omegas = omegas, 
                        p = p, 
                        beta0 = beta0, 
                        remove_isolates = TRUE)
```

### Fit JANE on network

``` r
res <- JANE::JANE(A = sim_data$A,
                  D = 2L,
                  K = 3L,
                  initialization = "GNN", 
                  model = "NDH",
                  case_control = FALSE,
                  DA_type = "none")
```

### Summarize and plot fit

``` r

# Summarize fit 
summary(res)

# plot network
plot(res)
```
