// helper_funs.h
#ifndef HELPER_FUNS_H
#define HELPER_FUNS_H

#include <RcppArmadillo.h>

// Inverse logit function
double logit_inv(double eta);

// Truncated Poisson density 
double trunc_poisson_density(double w, double mean, double log);

// Log-normal density 
double lognormal_density(double w, double precision, double mean, double log);

#endif // HELPER_FUNS_H