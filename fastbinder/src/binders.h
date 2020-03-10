/*
 * binders.h
 *
 *  Created on: Mar 6, 2020
 *      Author: toky
 */

#ifndef FASTBINDER_SRC_BINDERS_H_
#define FASTBINDER_SRC_BINDERS_H_

#include "math_utils.h"

typedef  arma::Mat<arma::uword> binder_mat_t ;

arma::uword  naive_binder (binder_mat_t CI , double Const_Binder = 0.5) ;
arma::uword  opt_binder (binder_mat_t CI , double Const_Binder = 0.5) ;
arma::uword  tiled_binder (binder_mat_t CI , double Const_Binder = 0.5) ;
arma::uword  parallel_binder (binder_mat_t CI , double Const_Binder = 0.5) ;
arma::uword  opencl_binder (binder_mat_t CI , double Const_Binder = 0.5) ;



#endif /* FASTBINDER_SRC_BINDERS_H_ */
