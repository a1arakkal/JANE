# Just Another Latent Space Network Clustering Algorithm (JANE)

# Overview

JANE is a package for fitting latent space network clustering models using an EM algorithm. 

# How to install

JANE can be installed via the `install_github()` function from the devtools package.

```
devtools::install_github("a1arakkal/JANE")
```

# Basic Usage

## Simulate a network

```
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

## Fit JANE on network

```
res <- JANE::JANE(A = sim_data$A,
                  D = 2L,
                  K = 3L,
                  initialization = "GNN", 
                  model = "NDH",
                  case_control = FALSE,
                  DA_type = "none")
```

## Summarise and plot fit

```
# Summarize fit 
summary(res)

# plot network
plot(res)
```

