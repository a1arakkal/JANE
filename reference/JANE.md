# Fit JANE

Fit a latent space cluster model, with or without noise edges, using an
EM algorithm.

## Usage

``` r
JANE(
  A,
  D = 2,
  K = 2,
  family = "bernoulli",
  noise_weights = FALSE,
  guess_noise_weights = NULL,
  model,
  initialization = "GNN",
  case_control = FALSE,
  DA_type = "none",
  seed = NULL,
  control = list()
)
```

## Arguments

- A:

  A square matrix or sparse matrix of class 'Matrix' (from the Matrix
  package) representing the adjacency matrix of the network of interest.

- D:

  Integer (scalar or vector) specifying the dimension of the latent
  space (default is 2).

- K:

  Integer (scalar or vector) specifying the number of clusters to
  consider (default is 2).

- family:

  A character string specifying the distribution of the edge weights.

  - 'bernoulli': Expects an **unweighted** network; utilizes a Bernoulli
    distribution with a logit link (default)

  - 'lognormal': Expects a **weighted** network with positive, non-zero,
    continuous edge weights; utilizes a log-normal distribution with an
    identity link

  - 'poisson': Expects a **weighted** network with edge weights
    representing non-zero counts; utilizes a zero-truncated Poisson
    distribution with a log link

- noise_weights:

  A logical; if `TRUE` then a Hurdle model is used to account for noise
  weights, if `FALSE` simply utilizes the supplied network (converted to
  an unweighted binary network if a weighted network is supplied, i.e.,
  (A \> 0.0)\*1.0) and fits a latent space cluster model (default is
  `FALSE`).

- guess_noise_weights:

  Only applicable if `noise_weights = TRUE`. A numeric value specifying
  the best guess for the mean of the noise weight distribution for
  `family %in% c('lognormal', 'poisson')` (mean is on the log-scale for
  lognormal) **OR** proportion (i.e. in (0,1)) of all edges that are
  noise edges for `family = 'bernoulli'`. If `NULL` (i.e., default) and
  `noise_weights = TRUE` then the 1st percentile of the non-zero weights
  will be used for `family %in% c('lognormal', 'poisson')` and 1% will
  be used for `family = 'bernoulli'`.

- model:

  A character string specifying the model to fit:

  - 'NDH': **undirected** network with no degree heterogeneity (or
    connection strength heterogeneity if working with weighted network)

  - 'RS': **undirected** network with degree heterogeneity (and
    connection strength heterogeneity if working with weighted network)

  - 'RSR': **directed** network with degree heterogeneity (and
    connection strength heterogeneity if working with weighted network)

- initialization:

  A character string or a list to specify the initial values for the EM
  algorithm:

  - 'GNN': uses a type of graphical neural network approach to generate
    initial values (default)

  - 'random': uses random initial values

  - A user-supplied list of S3
    [`class`](https://rdrr.io/r/base/class.html) "`JANE.initial_values`"
    representing initial values for the EM algorithm. See
    [`specify_initial_values`](https://a1arakkal.github.io/JANE/reference/specify_initial_values.md)
    for details on how to specify initial values

- case_control:

  A logical; if `TRUE` then uses a case-control approximation approach
  (default is `FALSE`).

- DA_type:

  (Experimental) A character string to specify the type of deterministic
  annealing approach to use

  - 'none': does not employ a deterministic annealing approach (default)

  - 'cooling': (Experimental) employs a traditional deterministic
    annealing approach where temperature decreases

  - 'heating': (Experimental) employs a deterministic anti-annealing
    approach where temperature increases

  - 'hybrid': (Experimental) employs a combination of the 'cooling' and
    'heating' approach

- seed:

  (optional) An integer value to specify the seed for reproducibility.

- control:

  A list of control parameters. See 'Details'.

## Value

A list of S3 [`class`](https://rdrr.io/r/base/class.html) "`JANE`"
containing the following components:

- `input_params`:

  A list containing the input parameters for `IC_selection`,
  `case_control`, `DA_type`, `model`, `family`, and `noise_weights` used
  in the function call.

- `A`:

  The square sparse adjacency matrix of class 'dgCMatrix' used in
  fitting the latent space cluster model. This matrix can be different
  than the input A matrix as isolates are removed.

- `IC_out`:

  A matrix containing the relevant information criteria for all
  combinations of `K`, `D`, and `n_start` considered. The 'selected'
  column indicates the chosen optimal fit.

- `all_convergence_ind`:

  A matrix containing the convergence information (i.e., 1 = converged,
  0 = did not converge) and number of iterations for all combinations of
  `K`, `D`, `n_start`, and `beta_temperature` considered.

- `optimal_res`:

  A list containing the estimated parameters of interest based on the
  optimal fit selected. It is recommended to use
  [`summary()`](https://rdrr.io/r/base/summary.html) to extract the
  parameters of interest. See
  [`summary.JANE`](https://a1arakkal.github.io/JANE/reference/summary.JANE.md)
  for more details.

- `optimal_starting`:

  A list of S3 [`class`](https://rdrr.io/r/base/class.html)
  "`JANE.initial_values`" containing the starting parameters used in the
  EM algorithm that resulted in the optimal fit selected. It is
  recommended to use [`summary()`](https://rdrr.io/r/base/summary.html)
  to extract the parameters of interest. See
  [`summary.JANE`](https://a1arakkal.github.io/JANE/reference/summary.JANE.md)
  for more details.

## Details

Isolates are removed from the adjacency matrix A. If an unsymmetric
adjacency matrix A is supplied for `model %in% c('NDH', 'RS')` the user
will be asked if they would like to proceed with converting A to a
symmetric matrix (i.e., `A <- 1.0 * ( (A + t(A)) > 0.0 )`); only able to
do so if `family = 'bernoulli'`. Additionally, if a weighted network is
supplied and `noise_weights = FALSE`, then the network will be converted
to an unweighted binary network (i.e., (A \> 0.0)\*1.0) and a latent
space cluster model is fit.

**`control`:**

The `control` argument is a named list that the user can supply
containing the following components:

- `verbose`:

  A logical; if `TRUE` causes additional information to be printed out
  about the progress of the EM algorithm (default is `FALSE`).

- `max_its`:

  An integer specifying the maximum number of iterations for the EM
  algorithm (default is `1e3`).

- `min_its`:

  An integer specifying the minimum number of iterations for the EM
  algorithm (default is `10`).

- `priors`:

  A list of S3 [`class`](https://rdrr.io/r/base/class.html)
  "`JANE.priors`" representing prior hyperparameter specifications
  (default is `NULL`). See
  [`specify_priors`](https://a1arakkal.github.io/JANE/reference/specify_priors.md)
  for details on how to specify the hyperparameters.

- `n_interior_knots`:

  (only relevant for `model %in% c('RS', 'RSR')`) An integer specifying
  the number of interior knots used in fitting a natural cubic spline
  for degree heterogeneity (and connection strength heterogeneity if
  working with weighted network) models (default is `5`).

- `termination_rule`:

  A character string to specify the termination rule to determine the
  convergence of the EM algorithm:

  - `'prob_mat'`: uses change in the absolute difference in
    \\\hat{Z}^{U}\\ (i.e., the \\N \times K\\ cluster membership
    probability matrix) between subsequent iterations (default)

  - `'Q'`: uses change in the absolute difference in the objective
    function of the E-step evaluated using parameters from subsequent
    iterations

  - `'ARI'`: comparing the classifications between subsequent iterations
    using adjusted Rand index

  - `'NMI'`: comparing the classifications between subsequent iterations
    using normalized mutual information

  - `'CER'`: comparing the classifications between subsequent iterations
    using classification error rate

- `tolerance`:

  A numeric specifying the tolerance used for
  `termination_rule %in% c('Q', 'prob_mat')` (default is `1e-3`).

- `tolerance_ARI`:

  A numeric specifying the tolerance used for `termination_rule = 'ARI'`
  (default is `0.999`).

- `tolerance_NMI`:

  A numeric specifying the tolerance used for `termination_rule = 'NMI'`
  (default is `0.999`).

- `tolerance_CER`:

  A numeric specifying the tolerance used for `termination_rule = 'CER'`
  (default is `0.01`).

- `n_its_start_CA`:

  An integer specifying what iteration to start computing the change in
  cumulative averages (note: the change in the cumulative average of
  \\\hat{U}\\, the latent position matrix, is not tracked when
  `termination_rule = 'Q'`) (default is `20`).

- `tolerance_diff_CA`:

  A numeric specifying the tolerance used for the change in cumulative
  average of `termination_rule` metric and \\\hat{U}\\ (note: the change
  in the cumulative average of \\\hat{U}\\ is not tracked when
  `termination_rule = 'Q'`) (default is `1e-3`).

- `consecutive_diff_CA`:

  An integer specifying the tolerance for the number of consecutive
  instances where the change in cumulative average is less than
  `tolerance_diff_CA` (default is `5`).

- `quantile_diff`:

  A numeric in `[0,1]` specifying the quantile used in computing the
  change in the absolute difference of \\\hat{Z}^{U}\\ and \\\hat{U}\\
  between subsequent iterations (default is `1`, i.e., max).

- `beta_temp_schedule`:

  (Experimental) A numeric vector specifying the temperature schedule
  for deterministic annealing (default is `1`, i.e., deterministic
  annealing not utilized).

- `n_control`:

  An integer specifying the fixed number of controls (i.e., non-links)
  sampled for each actor; only relevant when `case_control = TRUE`
  (default is `100` when `case_control = TRUE` and `NULL` when
  `case_control = FALSE`).

- `n_start`:

  An integer specifying the maximum number of starts for the EM
  algorithm (default is `5`).

- `max_retry`:

  An integer specifying the maximum number of re-attempts if starting
  values cause issues with EM algorithm (default is `5`).

- `IC_selection`:

  A character string to specify the information criteria used to select
  the optimal fit based on the combinations of `K`, `D`, and `n_start`
  considered:

  - `'BIC_model'`: BIC computed from logistic regression or Hurdle model
    component

  - `'BIC_mbc'`: BIC computed from model based clustering component

  - `'ICL_mbc'`: ICL computed from model based clustering component

  - `'Total_BIC'`: sum of `'BIC_model'` and `'BIC_mbc'`

  - `'Total_ICL'`: sum of `'BIC_model'` and `'ICL_mbc'` (default)

- `sd_random_U_GNN`:

  (only relevant when `initialization = 'GNN'`) A positive numeric value
  specifying the standard deviation for the random draws from a normal
  distribution to initialize \\U\\ (default is `1`).

- `max_retry_GNN`:

  (only relevant when `initialization = 'GNN'`) An integer specifying
  the maximum number of re-attempts for the `GNN` approach before
  switching to random starting values (default is `10`).

- `n_its_GNN`:

  (only relevant when `initialization = 'GNN'`) An integer specifying
  the maximum number of iterations for the `GNN` approach (default is
  `10`).

- `downsampling_GNN`:

  (only relevant when `initialization = 'GNN'`) A logical; if `TRUE`
  employs downsampling s.t. the number of links and non-links are
  balanced for the `GNN` approach (default is `TRUE`).

**Running `JANE` in parallel:**

`JANE` integrates the future and future.apply packages to fit the
various combinations of `K`, `D`, and `n_start` in parallel. The
'Examples' section below provides an example of how to run `JANE` in
parallel. See
[`plan`](https://future.futureverse.org/reference/plan.html) and
[`future.apply`](https://future.apply.futureverse.org/reference/future.apply.html)
for more details.

**Choosing the number of clusters:**

`JANE` allows for the following model selection criteria to choose the
number of clusters (smaller values are favored):

- `'BIC_model'`: BIC computed from logistic regression or Hurdle model
  component

- `'BIC_mbc'`: BIC computed from model based clustering component

- `'ICL_mbc'`: ICL (Biernacki et al. (2000)) computed from model based
  clustering component

- `'Total_BIC'`: Sum of `'BIC_model'` and `'BIC_mbc'`, this is the model
  selection criterion proposed by Handcock et al. (2007)

- `'Total_ICL'`: (default) sum of `'BIC_model'` and `'ICL_mbc'`, this
  criterion is similar to `'Total_BIC'`, but uses ICL for the model
  based clustering component, which tends to favor more well-separated
  clusters.

Based on simulation studies, Biernacki et al. (2000) recommends that
when the interest in mixture modeling is cluster analysis, instead of
density estimation, the \\ICL\_{mbc}\\ criterion should be favored over
the \\BIC\_{mbc}\\ criterion, as the \\BIC\_{mbc}\\ criterion tends to
overestimate the number of clusters. The \\BIC\_{mbc}\\ criterion is
designed to choose the number of components in a mixture model rather
than the number of clusters.

*Warning*: It is not certain whether it is appropriate to use the model
selection criterion above to select `D`.

## References

Biernacki, C., Celeux, G., Govaert, G., 2000. Assessing a mixture model
for clustering with the integrated completed likelihood. IEEE
Transactions on Pattern Analysis and Machine Intelligence 22, 719–725.

Handcock, M.S., Raftery, A.E., Tantrum, J.M., 2007. Model-based
clustering for social networks. Journal of the Royal Statistical Society
Series A: Statistics in Society 170, 301–354.

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

# Run JANE on simulated data - consider multiple D and K
res <- JANE::JANE(A = sim_data$A,
                  D = 2:5,
                  K = 2:10,
                  initialization = "GNN", 
                  model = "NDH",
                  case_control = FALSE,
                  DA_type = "none")
                  
# Run JANE on simulated data - parallel with 5 cores (NOT RUN)
# future::plan(future::multisession, workers = 5)
# res <- JANE::JANE(A = sim_data$A,
#                   D = 2L,
#                   K = 3L,
#                   initialization = "GNN", 
#                   model = "NDH",
#                   case_control = FALSE,
#                   DA_type = "none")
# future::plan(future::sequential)

# Run JANE on simulated data - case/control approach with 20 controls sampled for each actor
res <- JANE::JANE(A = sim_data$A,
                  D = 2L,
                  K = 3L,
                  initialization = "GNN", 
                  model = "NDH",
                  case_control = TRUE,
                  DA_type = "none",
                  control = list(n_control = 20))
                   
# Reproducibility
res1 <- JANE::JANE(A = sim_data$A,
                   D = 2L,
                   K = 3L,
                   initialization = "GNN", 
                   seed = 1234,
                   model = "NDH",
                   case_control = FALSE,
                   DA_type = "none")

res2 <- JANE::JANE(A = sim_data$A,
                   D = 2L,
                   K = 3L,
                   initialization = "GNN", 
                   seed = 1234,
                   model = "NDH",
                   case_control = FALSE,
                   DA_type = "none")  

## Check if results match
all.equal(res1, res2)    

# Another reproducibility example where the seed was not set. 
# It is possible to replicate the results using the starting values due to 
# the nature of EM algorithms
res3 <- JANE::JANE(A = sim_data$A,
                   D = 2L,
                   K = 3L,
                   initialization = "GNN", 
                   model = "NDH",
                   case_control = FALSE,
                   DA_type = "none")
## Extract starting values                    
start_vals <- res3$optimal_start

## Run JANE using extracted starting values, no need to specify D and K 
## below as function will determine those values from start_vals
res4 <- JANE::JANE(A = sim_data$A,
                   initialization = start_vals, 
                   model = "NDH",
                   case_control = FALSE,
                   DA_type = "none")
                   
## Check if optimal_res are identical
all.equal(res3$optimal_res, res4$optimal_res)                   
# }                            
```
