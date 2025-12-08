# JANE User Guide

## Description

**JANE** is an R package for fitting latent space network cluster models
using an expectation-maximization (EM) algorithm. It enables flexible
modeling of unweighted or weighted network data, with or without noise
edges, and supports both directed and undirected networks, with or
without degree and strength heterogeneity. Designed to efficiently
handle large networks, **JANE** allows users to explore latent
structure, identify actor-centric communities, and simulate networks
with customizable clustering and connectivity patterns.

Details on the methodology underlying the package can be found in
Arakkal and Sewell (2025).

## Features

- Fine-tuned algorithm for fast and accurate actor-centric clustering
  using latent space models
- Supports both weighted and unweighted networks, with or without noisy
  edges
- Allows for modeling of discrete and continuous edge weights
- Multiple initialization strategies for the EM algorithm, including
  random initialization, a novel graphical neural network
  initialization, and user-specified values
- Customizable prior hyperparameter specification
- Supports parallel implementation for faster inference using the
  **future** and **future.apply** packages (Bengtsson 2021)
- Supports case-control likelihood approximation for scalable inference
  (Raftery et al. 2012)
- Built-in visualization and summary tools

## Getting Started

### Installation

``` r
# Current release from CRAN
install.packages("JANE")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("a1arakkal/JANE")

# Development version with vignettes
devtools::install_github("a1arakkal/JANE", build_vignettes = TRUE)
```

### Documentation

Once installed, you can load the package and access its help
documentation using the following commands:

``` r
# Load JANE
library(JANE)

# Vignette
RShowDoc("JANE-User-Guide", package = "JANE")
```

## Fitting **JANE** on networks *without* noise edges

When the `noise_weights` argument in
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) is set to
`FALSE` (default), a standard latent space cluster model is fit to the
supplied *unweighted* network.

### Available models

#### `model = "NDH"`

Applicable for undirected networks and assumes no degree heterogeneity.
However, in real-world networks, it is rare to find no degree
heterogeneity, as most networks exhibit considerable variation in node
connectivity. Below is an example that fits an NDH model, specifying the
dimension of the latent space as 2 and the number of clusters as 3.

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 3,
     model = "NDH",
     noise_weights = FALSE)
```

#### `model = "RS"`

Applicable for undirected networks and adjusts for degree heterogeneity
through the inclusion of actor-specific sociality effects. Below is an
example that fits an RS model, where the number of clusters is specified
to be 3 and a range of dimensions between 2 and 5 for the latent space
is considered.

``` r
JANE(A = network_adjacency_matrix,
     D = 2:5,
     K = 3,
     model = "RS",
     noise_weights = FALSE)
```

#### `model = "RSR"`

Applicable for directed networks and adjusts for degree heterogeneity
through the inclusion of actor-specific sender and receiver effects.
Below is an example that fits an RSR model, where a range of cluster
numbers between 2 and 10 and a range of dimensions between 2 and 5 for
the latent space is considered.

``` r
JANE(A = network_adjacency_matrix,
     D = 2:5,
     K = 2:10,
     model = "RSR",
     noise_weights = FALSE)
```

## Fitting **JANE** on networks *with* noise edges

With [`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md),
*weighted* networks are particularly useful in settings where noise
edges are present. In such settings, when the `noise_weights` argument
in [`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) is set
to `TRUE`, a latent space hurdle model (LSHM) is fit. The LSHM leverages
information from both the propensity to form an edge and the observed
edge weights to probabilistically down-weight noisy edges, while
preserving edges that are structurally meaningful, subsequently
enhancing community detection.

If the supplied network is a weighted network, in the absence of noise
it can be shown that the latent positions, regression parameters
associated with the logistic regression model, finite mixture model
parameters, and actor-specific cluster membership probabilities can all
be estimated separately from the generalized linear model parameters for
the edge weights. **Thus, when `noise_weights = FALSE`,
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) will
simply dichotomize the supplied weighted network based on a threshold
value of 0 and fit a standard latent space cluster model**.

### Available models

You can fit the `"NDH"`, `"RS"`, or `"RSR"` models using any of the
supported weight distributions described below. When working with
weighted networks, the `"RS"` and `"RSR"` models account not only for
degree heterogeneity but also for strength heterogeneity.

#### `family = "poisson"`

Use this option for count-weighted networks. This setting models
observed edge weights using a zero-truncated Poisson (ZTP) distribution:

- Signal edge weights are modeled with a ZTP using a log link
- Noise edge weights are modeled with a ZTP using a fixed user-specified
  mean via the `guess_noise_weights` argument

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 5,
     model = "RS",
     noise_weights = TRUE,
     family = "poisson",
     guess_noise_weights = 1L)
```

#### `family = "lognormal"`

Use this option when edge weights are positive, continuous values. This
setting models observed edge weights using a log-normal distribution:

- Signal edges are modeled using a log-normal distribution with an
  identity link
- Noise edge weights are modeled using a log-normal distribution with a
  user-specified log-scale mean via the `guess_noise_weights` argument

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 5,
     model = "RS",
     noise_weights = TRUE,
     family = "lognormal",
     guess_noise_weights = -3.5) # log-scale mean
```

#### `family = "bernoulli"` (default)

This setting is used for unweighted (binary) networks. In this case:

- Signal and noise edges share the same observed weight (i.e., 1), so
  only the edge propensity model is used to identify and down-weight
  likely noise edges
- `guess_noise_weights` is interpreted as the expected proportion of all
  edges that are likely to be noise

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 5,
     model = "RSR",
     noise_weights = TRUE,
     family = "bernoulli",
     guess_noise_weights = 0.2) # expected noise edge proportion
```

### `guess_noise_weights`

If guess_noise_weights is left as NULL (the default),
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) will
automatically set this value based on the `family` argument:

- For `family %in% c("lognormal", "poisson")`, the 1st percentile of the
  non-zero edge weights is used (log-transformed weights for
  `lognormal`)
- For `family = "bernoulli"`, a default value of 0.01 (i.e., 1%) is
  used, representing the assumed proportion of noisy edges

## Simulating networks

**JANE** includes a built-in function,
[`sim_A()`](https://a1arakkal.github.io/JANE/reference/sim_A.md), for
simulating networks with varying clustering and connectivity patterns.
The function returns a list, which includes the following two
components:

- `A`: A sparse adjacency matrix of class ‘dgCMatrix’ representing the
  “true” underlying unweighted network with no noise edges
- `W`: A sparse adjacency matrix of class ‘dgCMatrix’ representing the
  unweighted or weighted network, with or without noise edges

The relationship between `A` and `W` depends on the values specified for
`family` and `noise_weights_prob`:

- If `family = "bernoulli"` and `noise_weights_prob = 0`, then `A` and
  `W` are identical (i.e., no noise simulated)
- If `family = "bernoulli"` and `noise_weights_prob > 0`, then `W`
  contains added noise edges, while `A` represents the noise-free
  underlying network
- If `family %in% c("poisson", "lognormal")` and
  `noise_weights_prob = 0`, then `W` contains continuous or count-valued
  edge weights, but dichotomizing `W` at 0 will reproduce `A`
- If `family %in% c("poisson", "lognormal")` and
  `noise_weights_prob > 0`, then `W` is a noisy, weighted network with
  additional spurious edges, and `A` remains the true underlying
  structure

Below is an example to simulate an unweighted, undirected, noise-free
network with a 2-dimensional latent space and 3 clusters, assuming no
degree heterogeneity.

``` r
# Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
mus <- matrix(c(-1,-1,  # Mean vector 1
                 1,-1,  # Mean vector 2
                 1, 1), # Mean vector 3
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)

# Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
                  diag(rep(7,2)),  # Precision matrix 2
                  diag(rep(7,2))), # Precision matrix 3
                  dim = c(2,2,3))

# Simulate a network
sim_A(N = 100L, 
      model = "NDH",
      mus = mus, 
      omegas = omegas, 
      p = rep(1/3, 3), 
      params_LR = list(beta0 = 1.0),
      remove_isolates = TRUE)
```

The parameter `beta0` above can be tuned to achieve a desired network
density:

``` r
# Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
mus <- matrix(c(-1,-1,  # Mean vector 1
                 1,-1,  # Mean vector 2
                 1, 1), # Mean vector 3
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)

# Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
                  diag(rep(7,2)),  # Precision matrix 2
                  diag(rep(7,2))), # Precision matrix 3
                  dim = c(2,2,3))

desired_density <- 0.1 # Target network density
min_density <- desired_density * 0.99  # Lower bound for acceptable density
max_density <- desired_density * 1.01  # Upper bound for acceptable density
n_act <- 100L # Number of actors in the network

density <- Inf # Initialize density to enter while loop
beta0 <- 0.5 # Initial value for intercept parameter
n_while_loop <- 0 # Counter for outer loop iterations
max_its <- 100 # Maximum number of iterations
change_beta0 <- 0.1  # Amount to adjust beta0 by

# Adjust beta0 until simulated network has the approximate desired density
while(! (density >= min_density & density <= max_density) ){
  
  if(n_while_loop>max_its){
    break
  }
  
  n_retry_isolate <- 0
  retry_isolate <- T
  
  # Retry until a network with no isolates is generated (this while loop is optional)
  while(retry_isolate){
    
    sim_data <- sim_A(N = n_act, 
                      model = "NDH",
                      mus = mus, 
                      omegas = omegas, 
                      p = rep(1/3, 3), 
                      params_LR = list(beta0 = beta0),
                      remove_isolates = TRUE)
    
    n_retry_isolate <- n_retry_isolate + 1
    
    # Accept network if no isolates remain, or if retried more than 10 times at the same beta0
    if(nrow(sim_data$A) == n_act | n_retry_isolate>10){
      retry_isolate <- F
    }
    
  }
  
  # Compute network density
  density <- igraph::graph.density(
    igraph::graph_from_adjacency_matrix(sim_data$A, mode = "undirected")
    )

  # Adjust beta0 based on density feedback
  if (density > max_density)  {
    beta0 <- beta0 - change_beta0
  }
  
  if (density < min_density)  {
    beta0 <- beta0 + change_beta0
  }
  
  n_while_loop <- n_while_loop + 1
  
}    

A <- sim_data$A # Final simulated adjacency matrix
igraph::graph.density(igraph::graph_from_adjacency_matrix(A, mode = "undirected")) # Verify density
```

Below is an example that simulates a directed, weighted network with
noise, degree and strength heterogeneity, a 2-dimensional latent space,
and 3 clusters.

``` r
# Specify the 3 x 2 matrix containing the 1 x 2 mean vectors of the 3 bivariate normals
mus <- matrix(c(-1,-1,  # Mean vector 1
                 1,-1,  # Mean vector 2
                 1, 1), # Mean vector 3
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)

# Specify the 2 x 2 x 3 array containing the 2 x 2 precision matrices of the 3 bivariate normals
omegas <- array(c(diag(rep(7,2)),  # Precision matrix 1
                  diag(rep(7,2)),  # Precision matrix 2
                  diag(rep(7,2))), # Precision matrix 3
                  dim = c(2,2,3))

# Simulate a network
sim_A(N = 100L, 
      model = "RSR",
      family = "poisson",
      mus = mus, 
      omegas = omegas, 
      p = rep(1/3, 3), 
      params_LR = list(beta0 = 1),
      params_weights = list(beta0 = 2),
      noise_weights_prob = 0.1,
      mean_noise_weights = 1,
      remove_isolates = TRUE)
```

## Model selection criteria for choosing the number of clusters

[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) allows
for the following model selection criteria to choose the number of
clusters and the best initialization of the EM algorithm (smaller values
are favored):

- `"Total_BIC"` (i.e., BIC): Sum of $BIC_{model}$ and $BIC_{mbc}$, where
  $BIC_{model}$ is the BIC computed from logistic regression or Hurdle
  model component and $BIC_{mbc}$ is the BIC computed from model based
  clustering component. This is the model selection criterion proposed
  by Handcock, Raftery, and Tantrum (2007)
- `"Total_ICL"` (i.e., BICL): (default) sum of $BIC_{model}$ and
  $ICL_{mbc}$, this criterion is similar to `"Total_BIC"`, but uses the
  integrated completed likelihood (ICL) for the model based clustering
  component, which tends to favor more well-separated clusters

Based on simulation studies, Biernacki, Celeux, and Govaert (2000)
recommends that when the interest in mixture modeling is cluster
analysis, instead of density estimation, the $ICL_{mbc}$ criterion
should be favored over the $BIC_{mbc}$ criterion, as the $BIC_{mbc}$
criterion tends to overestimate the number of clusters. The $BIC_{mbc}$
criterion is designed to choose the number of components in a mixture
model rather than the number of clusters.

***Warning***: It is not certain whether it is appropriate to use the
model selection criterion above to select the dimension of the latent
space $D$.

Below is an example that fits an RSR model, where a range of cluster
numbers between 2 and 10 is considered for a 2-dimensional latent space
and `"Total_BIC"` is used to select the optimal number of clusters and
best initialization of the EM algorithm.

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 2:10,
     model = "RSR",
     noise_weights = FALSE,
     control = list(IC_selection = "Total_BIC"))
```

Below is an example that fits an RSR model for a 2-dimensional latent
space with 3 clusters and `"Total_ICL"` is used to select the optimal
initialization of the EM algorithm from 10 unique starts. Note: the
number of starts for the EM algorithm is controlled through the
`control` argument by supplying a value for `n_start`.

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 3,
     model = "RSR",
     noise_weights = FALSE,
     control = list(IC_selection = "Total_ICL",
                    n_start = 10))
```

## Initialization of EM algorithm

[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) has three
initialization strategies to generate starting values for the EM
algorithm, controlled through the `initialization` argument:

- `"GNN"`: uses a type of graphical neural network approach to generate
  initial values (default). Details of the graphical neural network
  approach can be found in Arakkal and Sewell (2025)
- `"random"`: uses random initial values
- A user-supplied list of S3 class `"JANE.initial_values"` representing
  initial values for the EM algorithm. See
  [`?specify_initial_values`](https://a1arakkal.github.io/JANE/reference/specify_initial_values.md)
  for details on how to specify initial values

Below is an example where starting values are supplied to
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) using
[`specify_initial_values()`](https://a1arakkal.github.io/JANE/reference/specify_initial_values.md).

``` r
# Specify starting values
D <- 3
K <- 5
N <- nrow(sim_data$A)
n_interior_knots <- 5L
U <- matrix(stats::rnorm(N*D), nrow = N, ncol = D)
omegas <- stats::rWishart(n = K, df = D+1, Sigma = diag(D))
mus <- matrix(stats::rnorm(K*D), nrow = K, ncol = D)
p <- extraDistr::rdirichlet(n = 1, rep(3,K))[1,]
Z <-  extraDistr::rdirichlet(n = N, alpha = rep(1, K))
beta <- stats::rnorm(n = 1 + 2*(1 + n_interior_knots))

my_starting_values <- specify_initial_values(A = network_adjacency_matrix,
                                             D = D,
                                             K = K,
                                             model = "RSR",
                                             n_interior_knots = n_interior_knots,
                                             U = U,
                                             omegas = omegas, 
                                             mus = mus, 
                                             p = p, 
                                             Z = Z,
                                             beta = beta)         

# Run JANE using my_starting_values (no need to specify D and K as function will 
# determine those values from my_starting_values)
JANE(A = network_adjacency_matrix,
     initialization = my_starting_values,
     model = "RSR",
     noise_weights = FALSE)
```

## Specification of prior hyperparameters

### Prior distributions

The prior distributions are specified as follows:

#### Prior on ${\mathbf{μ}}_{k}$ and $\mathbf{\Omega}_{k}$

The same prior is used for $k = 1,\ldots,K$:

$$\mathbf{\Omega}_{k} \sim \text{Wishart}\left( c,\mathbf{G}^{- 1} \right)$$

$${\mathbf{μ}}_{k} \mid \mathbf{\Omega}_{k} \sim \text{MVN}\left( \mathbf{a},\left( b\mathbf{\Omega}_{k} \right)^{- 1} \right)$$

#### Prior on $\mathbf{p}$

For the current implementation, we require that all elements of the `nu`
vector be $\geq 1$ to prevent negative mixture weights for empty
clusters:

$$\mathbf{p} \sim \text{Dirichlet}\left( \nu_{1},\ldots,\nu_{K} \right)$$

#### Prior on ${\mathbf{β}}_{LR}$

$${\mathbf{β}}_{LR} \sim \text{MVN}\left( \mathbf{e},\mathbf{F}^{- 1} \right)$$

#### Prior on $q$

$$q \sim \text{Beta}(h,l)$$

#### *Zero-truncated Poisson*

##### Prior on ${\mathbf{β}}_{GLM}$

$${\mathbf{β}}_{GLM} \sim \text{MVN}\left( \mathbf{e}_{2},\mathbf{F}_{2}^{- 1} \right)$$

#### *Log-normal*

##### Prior on $\tau_{\text{weights}}^{2}$

$$\tau_{\text{weights}}^{2} \sim \text{Gamma}\left( \frac{m_{1}}{2},\frac{o_{1}}{2} \right)$$

##### Prior on ${\mathbf{β}}_{GLM}$

$${\mathbf{β}}_{GLM} \mid \tau_{\text{weights}}^{2} \sim \text{MVN}\left( \mathbf{e}_{2},\left( \tau_{\text{weights}}^{2}\mathbf{F}_{2} \right)^{- 1} \right)$$

##### Prior on $\tau_{\text{noise weights}}^{2}$

$$\tau_{\text{noise weights}}^{2} \sim \text{Gamma}\left( \frac{m_{2}}{2},\frac{o_{2}}{2} \right)$$

#### Default hyperparameters

- $\mathbf{a} = \mathbf{0}_{D}$
- $b = 1$
- $c = D + 1$
- $\mathbf{e} = \mathbf{0}_{dim{({\mathbf{β}}_{LR})}}$
- $\mathbf{e}_{2} = \mathbf{0}_{dim{({\mathbf{β}}_{GLM})}}$
- $\mathbf{F} = \text{Diag}\left( 1/100,\left( 1/2.5^{2} \right)\mathbf{1}_{dim{({\mathbf{β}}_{LR})} - 1} \right)$
- $\mathbf{F}_{2} = \text{Diag}\left( 1/100,\left( 1/2.5^{2} \right)\mathbf{1}_{dim{({\mathbf{β}}_{GLM})} - 1} \right)$
- $\mathbf{G} = \mathbf{I}_{D}$
- $h = l = 1$
- $m_{1} = m_{2} = o_{1} = o_{2} = 2$
- ${\mathbf{ν}} = 3\mathbf{1}_{K}$

where, $D$ is the dimension of the latent space, $K$ is the number of
clusters,
$$dim\left( {\mathbf{β}}_{LR} \right) = dim\left( {\mathbf{β}}_{GLM} \right) = \begin{cases}
1 & \text{"NDH”} \\
{\zeta + 2} & \text{"RS”} \\
{2\zeta + 3} & {\text{"RSR”},}
\end{cases}$$ and $\zeta$ is the number of interior knots used in
fitting a natural cubic spline for degree heterogeneity (and connection
strength heterogeneity if working with weighted network) models (i.e.,
`"RS"` and `"RSR"` only).

### Example of specifying hyperparameters for a single combination of $K$ and $D$

``` r
# Specify prior hyperparameters
D <- 3
K <- 5
n_interior_knots <- 5L
a <- rep(1, D)
b <- 3
c <- 4
G <- 10*diag(D)
nu <- rep(2, K)
e <- rep(0.5, 1 + (n_interior_knots + 1))
f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))

my_prior_hyperparameters <- specify_priors(D = D,
                                           K = K,
                                           model = "RS",
                                           n_interior_knots = n_interior_knots,
                                           a = a,
                                           b = b,
                                           c = c,
                                           G = G,
                                           nu = nu,
                                           e = e,
                                           f = f)
                                           
# Run JANE using supplied prior hyperparameters (no need to specify D and K 
# as function will determine those values from my_prior_hyperparameters)
JANE(A = network_adjacency_matrix,
     initialization = "GNN",
     model = "RS",
     noise_weights = FALSE,
     control = list(priors = my_prior_hyperparameters))
```

### Example of specifying hyperparameters for multiple combinations of $K$ and $D$

#### Nested list

``` r
# Specify first set of prior hyperparameters
D <- 3
K <- 5
n_interior_knots <- 5L
a <- rep(1, D)
b <- 3
c <- 4
G <- 10*diag(D)
nu <- rep(2, K)
e <- rep(0.5, 1 + (n_interior_knots + 1))
f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))

my_prior_hyperparameters_1 <- specify_priors(D = D,
                                             K = K,
                                             model = "RS",
                                             n_interior_knots = n_interior_knots,
                                             a = a,
                                             b = b,
                                             c = c,
                                             G = G,
                                             nu = nu,
                                             e = e,
                                             f = f)

# Specify second set of prior hyperparameters
D <- 2
K <- 3
n_interior_knots <- 5L
a <- rep(1, D)
b <- 3
c <- 4
G <- 10*diag(D)
nu <- rep(2, K)
e <- rep(0.5, 1 + (n_interior_knots + 1))
f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))

my_prior_hyperparameters_2 <- specify_priors(D = D,
                                             K = K,
                                             model = "RS",
                                             n_interior_knots = n_interior_knots,
                                             a = a,
                                             b = b,
                                             c = c,
                                             G = G,
                                             nu = nu,
                                             e = e,
                                             f = f)

# Create nested list
my_prior_hyperparameters <- list(my_prior_hyperparameters_1,
                                 my_prior_hyperparameters_2)

# Run JANE using supplied prior hyperparameters (no need to specify D and K 
# as function will determine those values from my_prior_hyperparameters)
JANE(A = network_adjacency_matrix,
     initialization = "GNN",
     model = "RS",
     noise_weights = FALSE,
     control = list(priors = my_prior_hyperparameters))
```

#### Unevaluated calls

Unevaluated calls using
[`quote()`](https://rdrr.io/r/base/substitute.html) can be supplied as
values for specific hyperparameters. This is particularly useful when
running [`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md)
for multiple combinations of $K$ and $D$.

``` r
# Specify prior hyperparameters as unevaluated calls
n_interior_knots <- 5L
e <- rep(0.5, 1 + (n_interior_knots + 1))
f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
my_prior_hyperparameters <- specify_priors(model = "RS",
                                           n_interior_knots = n_interior_knots,
                                           a = quote(rep(1, D)),
                                           b = 3,
                                           c = quote(D + 1),
                                           G = quote(10*diag(D)),
                                           nu = quote(rep(2, K)),
                                           e = e,
                                           f = f)
                                           
# Run JANE using supplied prior hyperparameters
JANE(A = network_adjacency_matrix,
     D = 2:5,
     K = 2:10,
     initialization = "GNN",
     model = "RS",
     noise_weights = FALSE,
     control = list(priors = my_prior_hyperparameters))
```

## Parallel implementation of `JANE()`

You can speed up fitting models for multiple combinations of cluster
number `K`, latent space dimension `D`, and number of initializations of
the EM algorithm `n_start`, by running
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) in
parallel. This leverages the `future` and `future.apply` packages to
distribute computation across cores.

Below is an example using 5 parallel workers to run
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) with
parallel backend enabled, followed by resetting to sequential
processing.

``` r
# Set up parallel plan with 5 workers (cores)
future::plan(future::multisession, workers = 5)

# Run JANE in parallel 
res_parallel <- JANE(A = network_adjacency_matrix,
                     D = 2:5,
                     K = 3:10,
                     initialization = "GNN",
                     model = "RSR")

# Reset to sequential processing
future::plan(future::sequential)
```

## Case-control implementation of `JANE()`

When working with large *sparse* networks, the case-control
approximation implemented in
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) can
reduce computational cost. This approach leverages approximation methods
developed by Raftery et al. (2012). With the case-control
implementation, the likelihood includes all edges (where most of the
information lies) but uses a Monte Carlo approximation of terms
involving the unconnected dyads of the network (i.e., non-edges). Thus,
the cost of evaluating the likelihood becomes linear, rather than
quadratic, in $N$.

Below is an example running
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) with the
case-control approach, sampling 20 controls (i.e., non-edges) per actor.

``` r
JANE(A = network_adjacency_matrix,
     D = 2,
     K = 3,
     initialization = "GNN",
     model = "RSR",
     case_control = TRUE,
     control = list(n_control = 20))
```

## Using S3 methods with **JANE** objects

Once you have fit a model using the
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) function,
you can inspect and analyze the results using several built-in S3
methods. These methods work with the `JANE` S3 class object returned
from a model fit.

For the examples in this section, we assume that `res` is the object
returned by
[`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md).

### `print()`

The [`print()`](https://rdrr.io/r/base/print.html) method gives a
summary of the optimal model fit and available components.

``` r
print(res)
```

This displays:

- The model type, number of clusters, and latent space dimension
- Sizes of each cluster, that are estimated using a hard cluster
  assignment rule defined as
  $\{ h|{\widehat{Z}}_{ih}^{U} = max_{k}{\widehat{Z}}_{ik}^{U}\}$, where
  ${\widehat{Z}}_{ik}^{U}$ is the estimated conditional probability that
  the $i^{th}$ actor belongs to the $k^{th}$ cluster. Note, here the
  optimal number of clusters selected may include empty clusters based
  on the hard cluster assignments
- Available components in the fitted object `res`

### `summary()`

Use the [`summary()`](https://rdrr.io/r/base/summary.html) method to
extract and display detailed information about the estimated model.

``` r
summary(res)
```

You can also compare estimated cluster assignments to known true
assignments using the `true_labels` argument, or inspect starting values
using `initial_values = TRUE`.

``` r
summary(res, true_labels = true_labels_vec)
summary(res, initial_values = TRUE)
```

Output includes:

- Coefficients from logistic or GLM components
- Latent positions, cluster means, and precision matrices
- Soft cluster memberships and hard cluster assignments
- Clustering uncertainty
- ARI, NMI, classification error rate (if `true_labels` is supplied)
- Model selection criteria
- Summary of noise edges if `noise_weights = TRUE`. This includes the
  estimated conditional probability that a specific edge is a noise or a
  non-noise edge, and their corresponding noise edge cluster assignment
  based on a hard clustering rule of
  $\{ h|{\widehat{Z}}_{eh}^{W} = max\left( {\widehat{Z}}_{e1}^{W},{\widehat{Z}}_{e2}^{W} \right)\}$
  for $e = 1\ldots,|E|$, where ${\widehat{Z}}_{e1}^{W}$ and
  ${\widehat{Z}}_{e2}^{W}$ are the estimated conditional probabilities
  that the $e^{th}$ edge is a non-noise and noise edge, respectively,
  and $|E|$ represents the total number of edges in the network (for
  undirected networks, only the upper diagonal edges are retained).
  Noise edge cluster labels are defined as: 1 = non-noise edge and 2 =
  noise edge

### `plot()`

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
provides visualizations depending on the `type` argument.

#### Latent space clustering (default `type = "lsnc"`)

``` r
plot(res)
```

- Plots latent positions colored by cluster
- Contours show the fitted Gaussian components when $D$ = 2. When $D$ \>
  2, pairwise plots of the dimensions are displayed instead

#### Misclassified actors (`type = "misclassified"`)

``` r
plot(res, type = "misclassified", true_labels = true_labels_vec)
```

- Highlights misclassified actors in black based on supplied true
  cluster assignments (i.e., `true_labels_vec` in this example)

#### Clustering uncertainty (`type = "uncertainty"`)

``` r
plot(res, type = "uncertainty")
```

- Colors nodes by uncertainty (i.e., 1 -
  $max_{k}{\widehat{Z}}_{ik}^{U}$)

#### EM trace plot (`type = "trace_plot"`)

``` r
plot(res, type = "trace_plot")
```

- Shows convergence metrics over EM iterations

#### Additional options

- `alpha_edge`, `alpha_node`, `zoom`, `swap_axes`, `rotation_angle`,
  `cluster_cols`, etc.
- `remove_noise_edges = TRUE` removes noise edges based on noise edge
  hard cluster assignment (only applicable if
  [`JANE()`](https://a1arakkal.github.io/JANE/reference/JANE.md) was run
  with `noise_weights = TRUE`)

## References

Arakkal, Alan T, and Daniel K Sewell. 2025. “JANE: Just Another Latent
Space NEtwork Clustering Algorithm.” *Computational Statistics & Data
Analysis* 211: 108228.

Bengtsson, Henrik. 2021. “A Unifying Framework for Parallel and
Distributed Processing in r Using Futures.” *The R Journal* 13 (2):
208–27. <https://doi.org/10.32614/RJ-2021-048>.

Biernacki, Christophe, Gilles Celeux, and Gérard Govaert. 2000.
“Assessing a Mixture Model for Clustering with the Integrated Completed
Likelihood.” *IEEE Transactions on Pattern Analysis and Machine
Intelligence* 22 (7): 719–25.

Handcock, Mark S, Adrian E Raftery, and Jeremy M Tantrum. 2007.
“Model-Based Clustering for Social Networks.” *Journal of the Royal
Statistical Society Series A: Statistics in Society* 170 (2): 301–54.

Raftery, Adrian E, Xiaoyue Niu, Peter D Hoff, and Ka Yee Yeung. 2012.
“Fast Inference for the Latent Space Network Model Using a Case-Control
Approximate Likelihood.” *Journal of Computational and Graphical
Statistics* 21 (4): 901–19.
