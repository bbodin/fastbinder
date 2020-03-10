/*
 * maths.hpp
 *
 *  Created on: Jun 14, 2019
 */

#ifndef FASTBINDER_MIXTURE_CPP_MATHS_HPP_
#define FASTBINDER_MIXTURE_CPP_MATHS_HPP_

#include "verbose.h"

#ifdef NO_RCPP
#include <armadillo>
#else
#ifdef HAS_RCPP
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#else
#error "Unsupported Compilation flags, Need NO_RCPP or HAS_RCPP"
#endif
#endif


#endif /* FASTBINDER_MIXTURE_CPP_MATHS_HPP_ */
