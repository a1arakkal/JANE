# Specify starting values for EM algorithm

A function that allows the user to specify starting values for the EM
algorithm in a structure accepted by
[`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md).

## Usage

``` r
specify_initial_values(
  A,
  D,
  K,
  model,
  family = "bernoulli",
  noise_weights = FALSE,
  n_interior_knots = NULL,
  U,
  omegas,
  mus,
  p,
  Z,
  beta,
  beta2,
  precision_weights,
  precision_noise_weights
)
```

## Arguments

- A:

  A square matrix or sparse matrix of class 'dgCMatrix' representing the
  adjacency matrix of the network of interest.

- D:

  An integer specifying the dimension of the latent positions.

- K:

  An integer specifying the total number of clusters.

- model:

  A character string specifying the model:

  - 'NDH': **undirected** network with no degree heterogeneity (or
    connection strength heterogeneity if working with weighted network)

  - 'RS': **undirected** network with degree heterogeneity (and
    connection strength heterogeneity if working with weighted network)

  - 'RSR': **directed** network with degree heterogeneity (and
    connection strength heterogeneity if working with weighted network)

- family:

  A character string specifying the distribution of the edge weights.

  - 'bernoulli': for **unweighted** networks; utilizes a Bernoulli
    distribution with a logit link (default)

  - 'lognormal': for **weighted** networks with positive, non-zero,
    continuous edge weights; utilizes a log-normal distribution with an
    identity link

  - 'poisson': for **weighted** networks with edge weights representing
    non-zero counts; utilizes a zero-truncated Poisson distribution with
    a log link

- noise_weights:

  A logical; if TRUE then a Hurdle model is used to account for noise
  weights, if FALSE simply utilizes the supplied network (converted to
  an unweighted binary network if a weighted network is supplied, i.e.,
  (A \> 0.0)\*1.0) and fits a latent space cluster model (default is
  FALSE).

- n_interior_knots:

  An integer specifying the number of interior knots used in fitting a
  natural cubic spline for degree heterogeneity (and connection strength
  heterogeneity if working with weighted network) models (i.e., 'RS' and
  'RSR' only; default is `NULL`).

- U:

  A numeric \\N \times D\\ matrix with rows specifying an actor's
  position in a \\D\\-dimensional social space.

- omegas:

  A numeric \\D \times D \times K\\ array specifying the precision
  matrices of the \\K\\ \\D\\-variate normal distributions for the
  latent positions.

- mus:

  A numeric \\K \times D\\ matrix specifying the mean vectors of the
  \\K\\ \\D\\-variate normal distributions for the latent positions.

- p:

  A numeric vector of length \\K\\ specifying the mixture weights of the
  finite multivariate normal mixture distribution for the latent
  positions.

- Z:

  A numeric \\N \times K\\ matrix with rows representing the conditional
  probability that an actor belongs to the cluster \\K = k\\ for \\k =
  1,\ldots,K\\.

- beta:

  A numeric vector specifying the regression coefficients for the
  logistic regression model. Specifically, a vector of length  
  `1 + (model =="RS")*(n_interior_knots + 1) +`  
  `(model =="RSR")*2*(n_interior_knots + 1)`.

- beta2:

  A numeric vector specifying the regression coefficients for the
  zero-truncated Poisson or log-normal GLM. Specifically, a vector of
  length  
  `1 + (model =="RS")*(n_interior_knots + 1) +`  
  `(model =="RSR")*2*(n_interior_knots + 1)`.  
  Only relevant when `noise_weights = TRUE & family != 'bernoulli'`.

- precision_weights:

  A positive numeric scalar specifying the precision (on the log scale)
  of the log-normal weight distribution. Only relevant when
  `noise_weights = TRUE & family = 'lognormal'`.

- precision_noise_weights:

  A positive numeric scalar specifying the precision (on the log scale)
  of the log-normal noise weight distribution. Only relevant when
  `noise_weights = TRUE & family = 'lognormal'`.

## Value

A list of S3 [`class`](https://rdrr.io/r/base/class.html)
"`JANE.initial_values`" representing starting values for the EM
algorithm, in a structure accepted by `JANE`.

## Details

To match [`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md),
this function will remove isolates from the adjacency matrix A and
determine the total number of actors after excluding isolates. If this
is not done, errors with respect to incorrect dimensions in the starting
values will be generated when executing
[`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md).

Similarly to match
[`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md), if an
unsymmetric adjacency matrix A is supplied for
`model %in% c('NDH', 'RS')` the user will be asked if they would like to
proceed with converting A to a symmetric matrix (i.e.,
`A <- 1.0 * ( (A + t(A)) > 0.0 )`). Additionally, if a weighted network
is supplied and `noise_weights = FALSE`, then the network will be
converted to an unweighted binary network (i.e., (A \> 0.0)\*1.0).

## Examples

``` r
# \donttest{
# Simulate network
mus <- matrix(c(-1,-1,1,-1,1,1), 
              nrow = 3,
              ncol = 2, 
              byrow = TRUE)
omegas <- array(c(diag(rep(7,2)),
                  diag(rep(7,2)), 
                  diag(rep(7,2))), 
                dim = c(2,2,3))
p <- rep(1/3, 3)
beta0 <- -1
sim_data <- JANE::sim_A(N = 100L, 
                        model = "RSR",
                        mus = mus, 
                        omegas = omegas, 
                        p = p, 
                        params_LR = list(beta0 = beta0),
                        remove_isolates = TRUE)

# Specify starting values
D <- 3L
K <- 5L
N <- nrow(sim_data$A)
n_interior_knots <- 5L

U <- matrix(stats::rnorm(N*D), nrow = N, ncol = D)
omegas <- stats::rWishart(n = K, df = D+1, Sigma = diag(D))
mus <- matrix(stats::rnorm(K*D), nrow = K, ncol = D)
p <- extraDistr::rdirichlet(n = 1, rep(3,K))[1,]
Z <-  extraDistr::rdirichlet(n = N, alpha = rep(1, K))
beta <- stats::rnorm(n = 1 + 2*(1 + n_interior_knots))

my_starting_values <- JANE::specify_initial_values(A = sim_data$A,
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
res <- JANE::JANE(A = sim_data$A,
                  initialization = my_starting_values,
                  model = "RSR")
# }
```
