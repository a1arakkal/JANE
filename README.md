
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JANE

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/JANE)](https://CRAN.R-project.org/package=JANE)
[![R-CMD-check](https://github.com/a1arakkal/JANE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a1arakkal/JANE/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/a1arakkal/JANE/branch/master/graph/badge.svg)](https://app.codecov.io/gh/a1arakkal/JANE?branch=master)

<!-- badges: end -->

**JANE** is an R package for fitting latent space network cluster models using an expectation-maximization (EM) algorithm. It enables flexible modeling of unweighted or weighted network data, with or without noise edges, and supports both directed and undirected networks, with or without degree and strength heterogeneity. Designed to efficiently handle large networks, **JANE** allows users to explore latent structure, identify actor-centric communities, and simulate networks with customizable clustering and connectivity patterns.

Details on the methodology underlying the package can be found [here](https://doi.org/10.1016/j.csda.2025.108228).

## Installation

``` r
# Current release from CRAN
install.packages("JANE")

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

# Plot network
plot(res)
```
