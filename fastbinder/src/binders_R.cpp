/*
 * binder.cpp
 *
 *  Created on: Mar 5, 2020
 *      Author: toky
 */


#include "binders.h"


// [[Rcpp::export]]
long IAM_binder_naive (Rcpp::IntegerMatrix  CI, int version = 0,   double Const_Binder = 0.5 ){

	arma::uword Binder_ind = naive_binder ( Rcpp::as<binder_mat_t>(CI) , Const_Binder);

	VERBOSE_DEBUG(" Shift by one for R index space ");
	Binder_ind++;


	return Binder_ind ;
}

// [[Rcpp::export]]
long IAM_binder_opt (Rcpp::IntegerMatrix  CI, int version = 0,   double Const_Binder = 0.5 ){

	arma::uword Binder_ind = opt_binder ( Rcpp::as<binder_mat_t>(CI) , Const_Binder);

	VERBOSE_DEBUG(" Shift by one for R index space ");
	Binder_ind++;


	return Binder_ind ;
}

// [[Rcpp::export]]
long IAM_binder_parallel (Rcpp::IntegerMatrix  CI, int version = 0,   double Const_Binder = 0.5 ){

	arma::uword Binder_ind = parallel_binder ( Rcpp::as<binder_mat_t>(CI) , Const_Binder);

	VERBOSE_DEBUG(" Shift by one for R index space ");
	Binder_ind++;


	return Binder_ind ;
}
