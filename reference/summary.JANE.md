# Summarizing JANE fits

S3 summary method for object of class "`JANE`".

## Usage

``` r
# S3 method for class 'JANE'
summary(object, true_labels = NULL, initial_values = FALSE, ...)
```

## Arguments

- object:

  An object of S3 [`class`](https://rdrr.io/r/base/class.html) "`JANE`",
  a result of a call to
  [`JANE`](https://a1arakkal.github.io/JANE/reference/JANE.md).

- true_labels:

  (optional) A numeric, character, or factor vector of known true
  cluster labels. Must have the same length as number of actors in the
  fitted network. Need to account for potential isolates removed
  (default is `NULL`).

- initial_values:

  A logical; if `TRUE` then summarize fit using the starting parameters
  used in the EM algorithm (default is `FALSE`, i.e., the results after
  the EM algorithm is run are summarized).

- ...:

  Unused.

## Value

A list of S3 [`class`](https://rdrr.io/r/base/class.html)
"`summary.JANE`" containing the following components (Note: \\N\\ is the
number of actors in the network, \\K\\ is the number of clusters, and
\\D\\ is the dimension of the latent space):

- `coefficients`:

  A list containing the estimated coefficients from the logistic
  regression model (i.e., 'beta_LR') and, if relevant, the estimated
  coefficients from the zero- truncated Poisson or log-normal GLM (i.e.,
  'beta_GLM').

- `U`:

  A numeric \\N \times D\\ matrix with rows containing an actor's
  estimated latent position in a \\D\\-dimensional social space.

- `p`:

  A numeric vector of length \\K\\ containing the estimated mixture
  weights of the finite multivariate normal mixture distribution for the
  latent positions.

- `mus`:

  A numeric \\K \times D\\ matrix containing the estimated mean vectors
  of the \\K\\ \\D\\-variate normal distributions for the latent
  positions.

- `omegas`:

  A numeric \\D \times D \times K\\ array containing the estimated
  precision matrices of the \\K\\ \\D\\-variate normal distributions for
  the latent positions.

- `Z_U`:

  A numeric \\N \times K\\ matrix with rows containing the estimated
  conditional probability that an actor belongs to the cluster \\K = k\\
  for \\k = 1,\ldots,K\\.

- `uncertainty`:

  A numeric vector of length \\N\\ containing the uncertainty of the
  \\i^{th}\\ actor's classification, derived as 1 - \\max_k
  \hat{Z}^{U}\_{ik}\\.

- `cluster_labels`:

  A numeric vector of length \\N\\ containing the cluster assignment of
  each actor based on a hard clustering rule of \\\\h \|
  \hat{Z}^{U}\_{ih} = max_k \hat{Z}^{U}\_{ik}\\\\.

- `Z_W`:

  A numeric \\\|E\| \times 6\\ matrix, with \\\|E\|\\ representing the
  total number of edges in the network (for undirected networks, only
  the upper diagonal edges are retained). The first two columns (i.e.,
  'i' and 'j') contains the specific indices of the edge between the
  \\i^{th}\\ and \\j^{th}\\ actors, the third column (i.e., 'weight')
  contains the specific edge weight, the fourth column (i.e.,
  'hat_zij1') contains the estimated conditional probability that the
  specific edge is a non-noise edge, the fifth column (i.e., 'hat_zij2')
  contains the estimated conditional probability that the specific edge
  is a noise edge, and the sixth column (i.e.,
  'noise_edge_cluster_labels') contains the noise-edge cluster
  assignment of each edge based on a hard clustering rule of \\\\h \|
  \hat{Z}^{W}\_{eh} = max(\hat{Z}^{W}\_{e1}, \hat{Z}^{W}\_{e2})\\\\ for
  \\e = 1,\ldots,\|E\|\\, where \\\hat{Z}^{W}\_{e1}\\ and
  \\\hat{Z}^{W}\_{e2}\\ are the estimated conditional probabilities that
  the \\e^{th}\\ edge is a non-noise and noise edge, respectively
  (labels defined as, 1: non-noise edge and 2: noise edge). Will be
  `NULL` if `noise_weights = FALSE` or `initial_values = TRUE`.

- `q_prob`:

  A numeric scalar representing the estimated proportion of non-edges in
  the "true" unobserved network that were converted to noise edges.

- `precision_weights`:

  A numeric scalar representing the estimated precision (on the log
  scale) of the log-normal weight distribution. Only relevant for
  `family = 'lognormal' & noise_weights = TRUE`.

- `precision_noise_weights`:

  A numeric scalar representing the estimated precision (on the log
  scale) of the log-normal noise weight distribution. Only relevant for
  `family = 'lognormal' & noise_weights = TRUE`.

- `IC`:

  Information criteria values of the optimal fit selected, including

  - `'BIC_model'`: BIC computed from logistic regression or Hurdle model
    component

  - `'BIC_mbc'`: BIC computed from model based clustering component

  - `'ICL_mbc'`: ICL computed from model based clustering component

  - `'Total_BIC'`: sum of `'BIC_model'` and `'BIC_mbc'`

  - `'Total_ICL'`: sum of `'BIC_model'` and `'ICL_mbc'`

- `input_params`:

  A list with the following components:

  - `model`: A character string containing the specific `model` used
    (i.e., 'NDH', 'RS', or 'RSR')

  - `family`: A character string containing the specific `family` used
    (i.e., 'bernoulli', 'poisson', or 'lognormal')

  - `noise_weights`: A logical; if `TRUE` then the approach utilizing a
    Hurdle model accounting for noise edges was utilized

  - `IC_selection`: A character string containing the specific
    information criteria used to select the optimal fit (i.e.,
    'BIC_model', 'BIC_mbc', 'ICL_mbc', 'Total_BIC', or 'Total_ICL')

  - `case_control`: A logical; if `TRUE` then the case/control approach
    was utilized

  - `DA_type`: A character string containing the specific deterministic
    annealing approach utilized (i.e., 'none', 'cooling', 'heating', or
    'hybrid')

  - `priors`: A list of the prior hyperparameters used. See
    [`specify_priors`](https://a1arakkal.github.io/JANE/reference/specify_priors.md)
    for definitions.

- `clustering_performance`:

  (only if `true_labels` is `!NULL`) A list with the following
  components:

  - `CER`: A list with two components: (i) `misclassified`: The indices
    of the misclassified actors in a minimum error mapping between the
    cluster labels and the known true cluster labels (i.e.,
    `true_labels`) and (ii) `errorRate`: The error rate corresponding to
    a minimum error mapping between the cluster labels and the known
    true cluster labels (see
    [`classError`](https://mclust-org.github.io/mclust/reference/classError.html)
    for details)

  - `ARI`: A numeric value containing the adjusted Rand index comparing
    the cluster labels and the known true cluster labels (see
    [`adjustedRandIndex`](https://mclust-org.github.io/mclust/reference/adjustedRandIndex.html)
    for details)

  - `NMI`: A numeric value containing the normalized mutual information
    comparing the cluster labels and the known true cluster labels (see
    [`NMI`](https://rdrr.io/pkg/aricode/man/NMI.html) for details)

  - `confusion_matrix`: A numeric table containing the confusion matrix
    comparing the cluster labels and the known true cluster labels.

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
                        model = "NDH",
                        mus = mus, 
                        omegas = omegas, 
                        p = p, 
                        params_LR = list(beta0 = beta0), 
                        remove_isolates = TRUE)
                        
# Run JANE on simulated data
res <- JANE::JANE(A = sim_data$A,
                  D = 2L,
                  K = 3L,
                  initialization = "GNN", 
                  model = "NDH",
                  case_control = FALSE,
                  DA_type = "none")
                  
# Summarize fit 
summary(res)

# Summarize fit and compare to true cluster labels
summary(res, true_labels = apply(sim_data$Z_U, 1, which.max))

# Summarize fit using starting values of EM algorithm
summary(res, initial_values = TRUE)
# }
```
