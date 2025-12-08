# Simulate unweighted or weighted networks, with or without noise edges, from latent space cluster models

Simulate an unweighted or weighted network, with or without noise edges,
from a \\D\\-dimensional latent space cluster model with \\K\\ clusters
and \\N\\ actors. The *squared* euclidean distance is used (i.e.,
\\dist(U_i,U_j)^2\\), where \\U_i\\ and \\U_j\\ are the respective
actor's positions in a \\D\\-dimensional social space.

## Usage

``` r
sim_A(
  N,
  mus,
  omegas,
  p,
  model = "NDH",
  family = "bernoulli",
  params_LR,
  params_weights = NULL,
  noise_weights_prob = 0,
  mean_noise_weights,
  precision_noise_weights,
  remove_isolates = TRUE
)
```

## Arguments

- N:

  An integer specifying the number of actors in the network.

- mus:

  A numeric \\K \times D\\ matrix specifying the mean vectors of the
  \\K\\ \\D\\-variate normal distributions for the latent positions.

- omegas:

  A numeric \\D \times D \times K\\ array specifying the precision
  matrices of the \\K\\ \\D\\-variate normal distributions for the
  latent positions.

- p:

  A numeric vector of length \\K\\ specifying the mixture weights of the
  finite multivariate normal mixture distribution for the latent
  positions.

- model:

  A character string specifying the type of model used to simulate the
  network:

  - 'NDH': generates an **undirected** network with no degree
    heterogeneity (or connection strength heterogeneity if working with
    weighted network)

  - 'RS': generates an **undirected** network with degree heterogeneity
    (and connection strength heterogeneity if working with weighted
    network), specifically by including actor specific random sociality
    effects

  - 'RSR': generates a **directed** network with degree heterogeneity
    (and connection strength heterogeneity if working with weighted
    network), specifically by including actor specific random sender and
    receiver effects

- family:

  A character string specifying the distribution of the edge weights.

  - 'bernoulli': generates an **unweighted** network from a latent space
    cluster model

  - 'lognormal': generates a **weighted** network by first simulating an
    unweighted network using a latent space cluster model, and then
    assigning edge weights based on a log-normal GLM utilizing an
    identity link

  - 'poisson': generates a **weighted** network by first simulating an
    unweighted network using a latent space cluster model, and then
    assigning edge weights based on a zero-truncated Poisson GLM
    utilizing a log link

- params_LR:

  A list containing the parameters of the logistic regression model to
  simulate the unweighted network, including:

  - 'beta0': a numeric value specifying the intercept parameter for the
    logistic regression model

  - 'precision_R_effects': precision parameters for random degree
    heterogeneity effects, specific to the logistic regression model:

    - 'NDH': does not apply, can leave as missing

    - 'RS': a numeric value specifying the precision parameter of the
      normal distribution of the random sociality effect, if missing
      will generate from a gamma(shape = 1, rate = 1)

    - 'RSR': a numeric matrix specifying the precision matrix of the
      multivariate normal distribution of the random sender and receiver
      effects, if missing will generate from a Wishart(df = 3, Sigma =
      \\I_2\\)

- params_weights:

  Only relevant when `family %in% c('lognormal', 'poisson')`. A list
  containing the parameters of the GLMs for the edge weights, including:

  - 'beta0': a numeric value specifying the intercept parameter for the
    zero-truncated Poisson or log-normal GLM

  - 'precision_R_effects': precision parameters for random connection
    strength heterogeneity effects, specific to the zero-truncated
    Poisson or log-normal GLM:

    - 'NDH': does not apply, can leave as missing

    - 'RS': a numeric value specifying the precision parameter of the
      normal distribution of the random sociality effect, if missing
      will generate from a gamma(shape = 1, rate = 1)

    - 'RSR': a numeric matrix specifying the precision matrix of the
      multivariate normal distribution of the random sender and receiver
      effects, if missing will generate from a Wishart(df = 3, Sigma =
      \\I_2\\)

  - 'precision_weights': a positive, non-zero, numeric representing the
    precision (on the log scale) of the log-normal weight distribution.
    Only relevant when `family = 'lognormal'`

- noise_weights_prob:

  A numeric in \[0,1\] representing the proportion of all edges in the
  simulated network that are noise edges (default is 0.0).

- mean_noise_weights:

  A numeric representing the mean of the noise weight distribution. Only
  relevant when `family %in% c('lognormal', 'poisson')` and
  `noise_weights_prob>0.0`. For family = 'poisson' value has to be \>
  0.0, for family = "lognormal" the mean is on the log scale.

- precision_noise_weights:

  A positive, non-zero, numeric representing the precision (on the log
  scale) of the log-normal noise weight distribution. Only relevant when
  `family = 'lognormal'` and `noise_weights_prob>0.0`.

- remove_isolates:

  A logical; if `TRUE` then isolates from the network are removed
  (default is `TRUE`).

## Value

A list containing the following components:

- `A`:

  A sparse adjacency matrix of class 'dgCMatrix' representing the "true"
  underlying unweighted network with no noise edges.

- `W`:

  A sparse adjacency matrix of class 'dgCMatrix' representing the
  unweighted or weighted network, with or without noise. Note, if
  `family = 'bernoulli'` and `noise_weights_prob = 0`, then `A = W`.

- `q_prob`:

  A numeric scalar representing the proportion of non-edges in the
  "true" underlying network converted to noise edges. See 'Details' for
  how this value is computed.

- `Z_U`:

  A numeric \\N \times K\\ cluster assignment matrix with rows
  representing the cluster an actor belongs to (i.e., indicated by a
  value of 1.0).

- `Z_W`:

  A numeric \\\|E\| \times 4\\ edge weight cluster assignment matrix,
  with \\\|E\|\\ representing the total number of edges in the network
  (for undirected networks, only the upper diagonal edges are retained).
  The first two columns (i.e., 'i' and 'j') contains the specific
  indices of the edge between the \\i^{th}\\ and \\j^{th}\\ actors, the
  third column (i.e., 'weight') contains the specific edge weight, and
  the fourth column (i.e., 'Z_W') represents a noise-cluster label,
  where 1 denotes a non-noise edge and 2 denotes a noise edge. Will be
  `NULL` if `noise_weights_prob = 0`.

- `U`:

  A numeric \\N \times D\\ matrix with rows representing an actor's
  position in a \\D\\-dimensional social space.

- `mus`:

  The inputted numeric \\K \times D\\ `mus` matrix.

- `omegas`:

  The inputted numeric \\D \times D \times K\\ `omegas` array.

- `p`:

  The inputted numeric vector `p` of length \\K\\.

- `noise_weights_prob`:

  The inputted numeric scalar `noise_weights_prob`.

- `mean_noise_weights`:

  The inputted numeric scalar `mean_noise_weights`. Will be `NULL` if
  `noise_weights_prob = 0`.

- `precision_noise_weights`:

  The inputted numeric scalar `precision_noise_weights`. Will be `NULL`
  if `noise_weights_prob = 0`.

- `model`:

  The inputted `model` character string.

- `family`:

  The inputted `family` character string.

- `params_LR`:

  The inputted `params_LR` list. If `model != "NDH"`, will have an
  additional element "RE" containing a numeric \\N \times 1\\ matrix
  representing the actor specific random sociality effect (i.e., s) OR a
  \\N \times 2\\ matrix representing the actor specific random sender
  and receiver effects (i.e., s and r, respectively).

- `params_weights`:

  The inputted `params_weights` list. If `model != "NDH"`, will have an
  additional element "RE" containing a numeric \\N \times 1\\ matrix
  representing the actor specific random sociality effect (i.e., s) OR a
  \\N \times 2\\ matrix representing the actor specific random sender
  and receiver effects (i.e., s and r, respectively).

## Details

The returned scalar `q_prob` represents the proportion of non-edges in
the simulated network to be converted to noise edges, computed as
\\\frac{p\_{noise} \times D\_{A}}{(1-D\_{A}) \times (1-p\_{noise})}\\,
where \\D\_{A}\\ is the density of the simulated network without noise
and \\p\_{noise}\\ is the inputted `noise_weights_prob`.

## Examples

``` r
# \donttest{

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

# Simulate an undirected, unweighted network, with no noise and no degree heterogeneity
JANE::sim_A(N = 100L, 
            model = "NDH",
            mus = mus, 
            omegas = omegas, 
            p = p, 
            params_LR = list(beta0 = beta0),
            remove_isolates = TRUE)

# Simulate a directed, weighted network, with degree and strength heterogeneity but no noise
JANE::sim_A(N = 100L, 
            model = "RSR",
            family = "lognormal",
            mus = mus, 
            omegas = omegas, 
            p = p, 
            params_LR = list(beta0 = beta0),
            params_weights = list(beta0 = 2,
                                  precision_weights = 1),
            remove_isolates = TRUE)

# Simulate an undirected, weighted network, with noise and degree and strength heterogeneity
JANE::sim_A(N = 100L, 
            model = "RS",
            family = "poisson",
            mus = mus, 
            omegas = omegas, 
            p = p, 
            params_LR = list(beta0 = beta0),
            params_weights = list(beta0 = 2),
            noise_weights_prob = 0.1,
            mean_noise_weights = 1,
            remove_isolates = TRUE)
# }
```
