
# JANE News

## Version Updates

## 0.1.0
* Initial public release

## 0.1.1
* Fix plotting scripts to include immediate call to on.exit to reset user settings
* Add testthat functionality for automatic testing of package functions

## 0.1.2
* Add import Rcpp and future to NAMESPACE

## 0.2.0
* Improve plotting function:
  * Add options for user-specified titles (i.e., main, xlab, and ylab)
  * Add option for user-specified colors for clusters 
  
## 0.2.1
* summary.JANE now returns actor-specific uncertainty  

## 1.0.0
* JANE is now able to handle weighted networks, specifically for applications that involve the presence of noise edges in the network

## 1.1.0
* Add importFrom progress progress_bar to NAMESPACE to address note in checks
* Add option to plotting function to remove noise edges if JANE was run with noise_weights = TRUE
* Add citation information

## 2.0.0
* Update prior specification functionality
* Make new S3 class for priors and initial values
* Improve documentation - add information about connection strength heterogeneity, fix priors in 'details' of specify_priors, and fix typos
* Add vignette "JANE User Guide"

## 2.0.1
* Add website with vignette to "URL" section of DESCRIPTION
* Fix bug in JANE when class is added to NULL optimal_starting
* Add more comprehensive tests
* Replace some helper functions with more numerically stable variants
* Use .registration = T for useDynLib
* @useDynLib only included once in zzz.R. Removed all other instances in other R scripts
* Create header file for helper functions
* Create helper functions for inverse logit, density for ZTP and log-normal, and solve_only_sympd (checks for sympd before running arma::solve)
* Update scripts to work with helper functions
* Fix deprecated methods::as(A, "dgCMatrix") in JANE.R
* Update specify_priors and check_priors to check user input, also remove redundant code
* Add dimension checks for user supplied mus, omegas, and p in sim_A
