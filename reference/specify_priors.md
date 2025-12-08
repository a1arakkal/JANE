# Specify prior hyperparameters for EM algorithm

A function that allows the user to specify the prior hyperparameters for
the EM algorithm in a structure accepted by
[`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md).

## Usage

``` r
specify_priors(
  D = 2,
  K = 2,
  model,
  family = "bernoulli",
  noise_weights = FALSE,
  n_interior_knots = NULL,
  a,
  b,
  c,
  G,
  nu,
  e,
  f,
  h,
  l,
  e_2,
  f_2,
  m_1,
  o_1,
  m_2,
  o_2
)
```

## Arguments

- D:

  An integer specifying the dimension of the latent positions (default
  is 2).

- K:

  An integer specifying the total number of clusters (default is 2).

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

- a:

  A numeric vector of length \\D\\ specifying the mean of the
  multivariate normal prior on \\\mu_k\\ for \\k = 1,\ldots,K\\, where
  \\\mu_k\\ represents the mean of the multivariate normal distribution
  for the latent positions of the \\k^{th}\\ cluster.

- b:

  A strictly positive numeric scalar specifying the scaling factor on
  the precision of the multivariate normal prior on \\\mu_k\\ for \\k =
  1,\ldots,K\\, where \\\mu_k\\ represents the mean of the multivariate
  normal distribution for the latent positions of the \\k^{th}\\
  cluster.

- c:

  A positive numeric scalar \\\> D\\ specifying the degrees of freedom
  of the Wishart prior on \\\Omega_k\\ for \\k = 1,\ldots,K\\, where
  \\\Omega_k\\ represents the precision of the multivariate normal
  distribution for the latent positions of the \\k^{th}\\ cluster.

- G:

  A numeric positive definite \\D \times D\\ matrix specifying the
  inverse of the scale matrix of the Wishart prior on \\\Omega_k\\ for
  \\k = 1,\ldots,K\\, where \\\Omega_k\\ represents the precision of the
  multivariate normal distribution for the latent positions of the
  \\k^{th}\\ cluster.

- nu:

  A positive numeric vector of length \\K\\ with elements \\\ge 1\\
  specifying the concentration parameters of the Dirichlet prior on
  \\p\\, where \\p\\ represents the mixture weights of the finite
  multivariate normal mixture distribution for the latent positions.

- e:

  A numeric vector of length
  `1 + (model =='RS')*(n_interior_knots + 1) + (model =='RSR')*2*(n_interior_knots + 1)`
  specifying the mean of the multivariate normal prior on
  \\\beta\_{LR}\\, where \\\beta\_{LR}\\ represents the coefficients of
  the logistic regression model.

- f:

  A numeric positive definite square matrix of dimension
  `1 + (model =='RS')*(n_interior_knots + 1) + (model =='RSR')*2*(n_interior_knots + 1)`
  specifying the precision of the multivariate normal prior on
  \\\beta\_{LR}\\, where \\\beta\_{LR}\\ represents the coefficients of
  the logistic regression model.

- h:

  A positive numeric scalar \\\ge 1\\ specifying the first shape
  parameter for the Beta prior on \\q\\, where \\q\\ is the proportion
  of non-edges in the "true" underlying network converted to noise
  edges. Only relevant when `noise_weights = TRUE`.

- l:

  A strictly positive numeric scalar specifying the second shape
  parameter for the Beta prior on \\q\\, where \\q\\ is the proportion
  of non-edges in the "true" underlying network converted to noise
  edges. Only relevant when `noise_weights = TRUE`.

- e_2:

  A numeric vector of length
  `1 + (model =='RS')*(n_interior_knots + 1) + (model =='RSR')*2*(n_interior_knots + 1)`
  specifying the mean of the multivariate normal prior on
  \\\beta\_{GLM}\\, where \\\beta\_{GLM}\\ represents the coefficients
  of the zero-truncated Poisson or log-normal GLM. Only relevant when
  `noise_weights = TRUE & family != 'bernoulli'`.

- f_2:

  A numeric positive definite square matrix of dimension
  `1 + (model =='RS')*(n_interior_knots + 1) + (model =='RSR')*2*(n_interior_knots + 1)`
  specifying the precision of the multivariate normal prior on
  \\\beta\_{GLM}\\, where \\\beta\_{GLM}\\ represents the coefficients
  of the zero-truncated Poisson or log-normal GLM. Only relevant when
  `noise_weights = TRUE & family != 'bernoulli'`.

- m_1:

  A positive numeric scalar \\\> 1\\ specifying the shape parameter for
  the Gamma prior on \\\tau^2\_{weights}\\, where \\\tau^2\_{weights}\\
  is the precision (on the log scale) of the log-normal weight
  distribution. Note, this value is scaled by 0.5, see 'Details'. Only
  relevant when `noise_weights = TRUE & family = 'lognormal'`.

- o_1:

  A positive numeric scalar \\\> 0\\ specifying the rate parameter for
  the Gamma prior on \\\tau^2\_{weights}\\, where \\\tau^2\_{weights}\\
  is the precision (on the log scale) of the log-normal weight
  distribution. Note, this value is scaled by 0.5, see 'Details'. Only
  relevant when `noise_weights = TRUE & family = 'lognormal'`.

- m_2:

  A positive numeric scalar \\\ge 2\\ specifying the shape parameter for
  the Gamma prior on \\\tau^2\_{noise \\ weights}\\, where
  \\\tau^2\_{noise \\ weights}\\ is the precision (on the log scale) of
  the log-normal noise weight distribution. Note, this value is scaled
  by 0.5, see 'Details'. Only relevant when
  `noise_weights = TRUE & family = 'lognormal'`.

- o_2:

  A positive numeric scalar \\\> 0\\ specifying the rate parameter for
  the Gamma prior on \\\tau^2\_{noise \\ weights}\\, where
  \\\tau^2\_{noise \\ weights}\\ is the precision (on the log scale) of
  the log-normal noise weight distribution. Note, this value is scaled
  by 0.5, see 'Details'. Only relevant when
  `noise_weights = TRUE & family = 'lognormal'`.

## Value

A list of S3 [`class`](https://rdrr.io/r/base/class.html)
"`JANE.priors`" representing prior hyperparameters for the EM algorithm,
in a structure accepted by `JANE`.

## Details

**Prior on \\\boldsymbol{\mu}\_k\\ and \\\boldsymbol{\Omega}\_k\\**
(note: the same prior is used for \\k = 1,\ldots,K\\) :

\$\$\boldsymbol{\Omega}\_k \sim Wishart(c, \boldsymbol{G}^{-1})\$\$
\$\$\boldsymbol{\mu}\_k \| \boldsymbol{\Omega}\_k \sim
MVN(\boldsymbol{a}, (b\boldsymbol{\Omega}\_k)^{-1})\$\$

**Prior on \\\boldsymbol{p}\\**:

For the current implementation we require that all elements of the `nu`
vector be \\\ge 1\\ to prevent against negative mixture weights for
empty clusters. \$\$\boldsymbol{p} \sim Dirichlet(\nu_1
,\ldots,\nu_K)\$\$

**Prior on \\\boldsymbol{\beta}\_{LR}\\**: \$\$\boldsymbol{\beta}\_{LR}
\sim MVN(\boldsymbol{e}, \boldsymbol{F}^{-1})\$\$

**Prior on \\q\\**: \$\$q \sim Beta(h, l)\$\$

***Zero-truncated Poisson***

**Prior on \\\boldsymbol{\beta}\_{GLM}\\**:
\$\$\boldsymbol{\beta}\_{GLM} \sim MVN(\boldsymbol{e}\_{2},
\boldsymbol{F}\_{2}^{-1})\$\$

***Log-normal***

**Prior on \\\tau^2\_{weights}\\**: \$\$\tau^2\_{weights} \sim
Gamma(\frac{m_1}{2}, \frac{o_1}{2})\$\$

**Prior on \\\boldsymbol{\beta}\_{GLM}\\**:
\$\$\boldsymbol{\beta}\_{GLM}\|\tau^2\_{weights} \sim
MVN(\boldsymbol{e}\_{2},
(\tau^2\_{weights}\boldsymbol{F}\_{2})^{-1})\$\$

**Prior on \\\tau^2\_{noise \\ weights}\\**: \$\$\tau^2\_{noise \\
weights} \sim Gamma(\frac{m_2}{2}, \frac{o_2}{2})\$\$

Unevaluated calls can be supplied as values for specific
hyperparameters. This is particularly useful when running `JANE` for
multiple combinations of `K` and `D`. See 'examples' section below for
implementation examples.

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
beta0 <- 1.0
sim_data <- JANE::sim_A(N = 100L, 
                        model = "RS",
                        mus = mus, 
                        omegas = omegas, 
                        p = p, 
                        params_LR = list(beta0 = beta0), 
                        remove_isolates = TRUE)
                        
                        
# Specify prior hyperparameters
D <- 3L
K <- 5L
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
                                           
# Run JANE on simulated data using supplied prior hyperparameters
res <- JANE::JANE(A = sim_data$A,
                  D = D,
                  K = K,
                  initialization = "GNN",
                  model = "RS",
                  case_control = FALSE,
                  DA_type = "none",
                  control = list(priors = my_prior_hyperparameters))

# Specify prior hyperparameters as unevaluated calls
n_interior_knots <- 5L
e <- rep(0.5, 1 + (n_interior_knots + 1))
f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))

my_prior_hyperparameters <- specify_priors(model = "RS",
                                           n_interior_knots = n_interior_knots,
                                           a = quote(rep(1, D)),
                                           b = b,
                                           c = quote(D + 1),
                                           G = quote(10*diag(D)),
                                           nu = quote(rep(2, K)),
                                           e = e,
                                           f = f)
                                           
# # Run JANE on simulated data using supplied prior hyperparameters (NOT RUN)
# future::plan(future::multisession, workers = 5)
# res <- JANE::JANE(A = sim_data$A,
#                    D = 2:5,
#                    K = 2:10,
#                    initialization = "GNN",
#                    model = "RS",
#                    case_control = FALSE,
#                    DA_type = "none",
#                    control = list(priors = my_prior_hyperparameters))
# future::plan(future::sequential)
                
                                                         
# }
```
