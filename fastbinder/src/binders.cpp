/*
 * binder.cpp
 *
 *  Created on: Mar 5, 2020
 *      Author: toky
 */

#include "binders.h"

#define ceild(n,d)  std::ceil(((double)(n))/((double)(d)))
#define floord(n,d) std::floor(((double)(n))/((double)(d)))

#define TILE_SIZE 64

/**
 * Naive implementation as a baseline
 */
arma::uword  naive_binder (arma::Mat<arma::uword> CI , double Const_Binder) {

	VERBOSE_DEBUG("Start seq_binder ");


	long N    = CI.n_cols;
	long iter = CI.n_rows;




	/**
	 * ## Compute similarity matrix pij ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){ pij <- pij + outer(c_out[[g]], c_out[[g]], "==")}
	 *
	 */


	VERBOSE_DEBUG("Init pij");

	arma::mat pij (N,N,arma::fill::zeros);

	for (int x = 0 ; x < N ; x++) {
		for (int y = 0 ; y < N ; y++) {
			pij(x,y) = 0;
			for (int d = 0 ; d < iter ; d++) {
				pij(x,y) += (CI(d,y) == CI(d,x)) ;
			}
	    }
	}


	VERBOSE_DEBUG("pij = pij/iter;");

	pij = pij/iter;

	/**
	 * ## Compute Binder_f              ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){
	 *                     cij <- outer(c_out[[g]], c_out[[g]], "==");
	 *                     aux <- (pij - Const_Binder) * as.matrix(cij);
	 *                     aux <-  aux[upper.tri(aux)];
	 *                     Binder_f[g] <- sum(aux);
	 *                  }
	 *
	 */

	VERBOSE_DEBUG(" Compute Binder_f    ");

	arma::vec Binder_f (iter,arma::fill::zeros);


	for (int cur_iter = 0 ; cur_iter < iter ; cur_iter++) {

		arma::mat cij (N,N,arma::fill::zeros);

		for (int x = 0 ; x < N ; x++) {
			for (int y = 0 ; y < N ; y++) {
				cij(x,y) = 0;
				cij(x,y) = (CI(cur_iter,y) == CI(cur_iter,x)) ;
			}
		}

		arma::mat aux1 = (pij - Const_Binder) % cij;
		arma::mat aux2 = arma::trimatu(aux1);
		Binder_f(cur_iter) = arma::accu(aux2);

	}


	/**
	 * ## Pick the maximum              ##
	 * ###################################
	 * Equivalent to:  	Binder_ind <- which.max(Binder_f)
	 *                  Binder_out <- c_out[Binder_ind]
	 *
	 */

	VERBOSE_DEBUG(" Return the maximum  ");
	arma::uword Binder_ind =  arma::index_max(Binder_f) ;


	return Binder_ind;

}



/**
 * Naive implementation as a baseline
 */
arma::uword  opt_binder (arma::Mat<arma::uword> CI , double Const_Binder ) {

	VERBOSE_DEBUG("Start seq_binder ");


	long N    = CI.n_cols;
	long iter = CI.n_rows;




	/**
	 * ## Compute similarity matrix pij ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){ pij <- pij + outer(c_out[[g]], c_out[[g]], "==")}
	 *
	 */


	VERBOSE_DEBUG("Init pij");

	arma::mat pij (N,N,arma::fill::zeros);

	for (int d = 0 ; d < iter ; d++) {
	for (int x = 0 ; x < N ; x++) {
		for (int y = x ; y < N ; y++) {
				pij(x,y) += (CI(d,y) == CI(d,x)) ;
			}
	    }
	}


	VERBOSE_DEBUG("pij = pij/iter;");

	pij = pij/iter;

	/**
	 * ## Compute Binder_f              ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){
	 *                     cij <- outer(c_out[[g]], c_out[[g]], "==");
	 *                     aux <- (pij - Const_Binder) * as.matrix(cij);
	 *                     aux <-  aux[upper.tri(aux)];
	 *                     Binder_f[g] <- sum(aux);
	 *                  }
	 *
	 */

	VERBOSE_DEBUG(" Compute Binder_f    ");

	arma::vec Binder_f (iter,arma::fill::zeros);

	arma::mat tmp = pij - Const_Binder;


		// double acc = 0.0;

	for (int cur_iter = 0 ; cur_iter < iter ; cur_iter++) {
	for (int x = 0 ; x < N ; x++) {
		for (int y = x ; y < N ; y++) {
				if (CI(cur_iter,y) == CI(cur_iter,x))
					Binder_f(cur_iter) += tmp(x,y)   ;
			}
		}
	}


	/**
	 * ## Pick the maximum              ##
	 * ###################################
	 * Equivalent to:  	Binder_ind <- which.max(Binder_f)
	 *                  Binder_out <- c_out[Binder_ind]
	 *
	 */

	VERBOSE_DEBUG(" Return the maximum  ");
	arma::uword Binder_ind =  arma::index_max(Binder_f) ;


	return Binder_ind;

}


/**
 * Parallel implementation as a baseline
 */
arma::uword  parallel_binder (arma::Mat<arma::uword> CI , double Const_Binder ) {

	VERBOSE_DEBUG("Start seq_binder ");


	int N    = CI.n_cols;
	int iter = CI.n_rows;




	/**
	 * ## Compute similarity matrix pij ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){ pij <- pij + outer(c_out[[g]], c_out[[g]], "==")}
	 *
	 */


	VERBOSE_DEBUG("Init pij");

	arma::mat pij (N,N,arma::fill::zeros);

	{
	  int t1, t2, t3, t4, t5, t6;
	 int lb, ub, lbp, ubp, lb2, ub2;
	 register int lbv, ubv;
	if ((N >= 1) && (iter >= 1)) {
	  lbp=0;
	  ubp=floord(N-1,TILE_SIZE);
	#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
	  for (t1=lbp;t1<=ubp;t1++) {
	    for (t2=t1;t2<=floord(N-1,TILE_SIZE);t2++) {
	      for (t3=0;t3<=floord(iter-1,TILE_SIZE);t3++) {
	        for (t4=TILE_SIZE*t1;t4<=std::min(N-1,TILE_SIZE*t1+(TILE_SIZE-1));t4++) {
	          for (t5=TILE_SIZE*t3;t5<=(std::min(iter-1,TILE_SIZE*t3+(TILE_SIZE-1)))-7;t5+=8) {
	            lbv=std::max(TILE_SIZE*t2,t4);
	            ubv=std::min(N-1,TILE_SIZE*t2+(TILE_SIZE-1));
	#pragma ivdep
	#pragma vector always
	            for (t6=lbv;t6<=ubv;t6++) {
	            	pij(t4,t6) += (CI(t5 + 0,t6) == CI(t5 + 0,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 1,t6) == CI(t5 + 1,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 2,t6) == CI(t5 + 2,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 3,t6) == CI(t5 + 3,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 4,t6) == CI(t5 + 4,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 5,t6) == CI(t5 + 5,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 6,t6) == CI(t5 + 6,t4)) ;
	            	pij(t4,t6) += (CI(t5 + 7,t6) == CI(t5 + 7,t4)) ;
	            }
	          }
	          for (;t5<=std::min(iter-1,TILE_SIZE*t3+(TILE_SIZE-1));t5++) {
	            lbv=std::max(TILE_SIZE*t2,t4);
	            ubv=std::min(N-1,TILE_SIZE*t2+(TILE_SIZE-1));
	#pragma ivdep
	#pragma vector always
	            for (t6=lbv;t6<=ubv;t6++) {
	              pij(t4,t6) += (CI(t5,t6) == CI(t5,t4)) ;
	            }
	          }
	        }
	      }
	    }
	  }
	}

}


	VERBOSE_DEBUG("pij = pij/iter;");

	pij = pij/iter;

	/**
	 * ## Compute Binder_f              ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){
	 *                     cij <- outer(c_out[[g]], c_out[[g]], "==");
	 *                     aux <- (pij - Const_Binder) * as.matrix(cij);
	 *                     aux <-  aux[upper.tri(aux)];
	 *                     Binder_f[g] <- sum(aux);
	 *                  }
	 *
	 */

	VERBOSE_DEBUG(" Compute Binder_f    ");

	arma::vec Binder_f (iter,arma::fill::zeros);

	arma::mat tmp = pij - Const_Binder;
	{
	  int t1, t2, t3, t4, t5;
	 int lb, ub, lbp, ubp, lb2, ub2;
	 register int lbv, ubv;
	if ((N >= 1) && (iter >= 1)) {
	  lbp=0;
	  ubp=floord(iter-1,TILE_SIZE);
	#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5)
	  for (t1=lbp;t1<=ubp;t1++) {
	    for (t2=0;t2<=floord(N-1,TILE_SIZE);t2++) {
	      for (t3=TILE_SIZE*t2;t3<=std::min(N-1,TILE_SIZE*t2+(TILE_SIZE-1));t3++) {
	        for (t4=t3;t4<=N-1-7;t4+=8) {
	          lbv=TILE_SIZE*t1;
	          ubv=std::min(iter-1,TILE_SIZE*t1+(TILE_SIZE-1));
	#pragma ivdep
	#pragma vector always
	          for (t5=lbv;t5<=ubv;t5++) {

	            if (CI(t5,t4 + 0) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 0)   ;
	            if (CI(t5,t4 + 1) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 1)   ;
	            if (CI(t5,t4 + 2) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 2)   ;
	            if (CI(t5,t4 + 3) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 3)   ;
	            if (CI(t5,t4 + 4) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 4)   ;
	            if (CI(t5,t4 + 5) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 5)   ;
	            if (CI(t5,t4 + 6) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 6)   ;
	            if (CI(t5,t4 + 7) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4 + 7)   ;
	          }
	        }
	        for (;t4<=N-1;t4++) {
	          lbv=TILE_SIZE*t1;
	          ubv=std::min(iter-1,TILE_SIZE*t1+(TILE_SIZE-1));
	#pragma ivdep
	#pragma vector always
	          for (t5=lbv;t5<=ubv;t5++) {
	            if (CI(t5,t4) == CI(t5,t3)) Binder_f(t5) += tmp(t3,t4)   ;
	          }
	        }
	      }
	    }
	  }
	}

}



	/**
	 * ## Pick the maximum              ##
	 * ###################################
	 * Equivalent to:  	Binder_ind <- which.max(Binder_f)
	 *                  Binder_out <- c_out[Binder_ind]
	 *
	 */

	VERBOSE_DEBUG(" Return the maximum  ");
	arma::uword Binder_ind =  arma::index_max(Binder_f) ;


	return Binder_ind;

}


/**
 * Tiled implementation as a baseline
 */
arma::uword  tiled_binder (arma::Mat<arma::uword> CI , double Const_Binder ) {

	VERBOSE_DEBUG("Start seq_binder ");


	int N    = CI.n_cols;
	int iter = CI.n_rows;




	/**
	 * ## Compute similarity matrix pij ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){ pij <- pij + outer(c_out[[g]], c_out[[g]], "==")}
	 *
	 */


	VERBOSE_DEBUG("Init pij");

	arma::mat pij (N,N,arma::fill::zeros);

	  int t1, t2, t3, t4, t5, t6;
	 register int lbv, ubv;
	if ((N >= 1) && (iter >= 1)) {
	  for (t1=0;t1<=floord(N-1,TILE_SIZE);t1++) {
	    for (t2=t1;t2<=floord(N-1,TILE_SIZE);t2++) {
	      for (t3=0;t3<=floord(iter-1,TILE_SIZE);t3++) {
	        for (t4=TILE_SIZE*t1;t4<=std::min(N-1,TILE_SIZE*t1+(TILE_SIZE-1));t4++) {
	          for (t5=TILE_SIZE*t3;t5<=std::min(iter-1,TILE_SIZE*t3+(TILE_SIZE-1));t5++) {
	            lbv=std::max(TILE_SIZE*t2,t4);
	            ubv=std::min(N-1,TILE_SIZE*t2+(TILE_SIZE-1));
	            for (t6=lbv;t6<=ubv;t6++) {
	              pij(t4,t6) += (CI(t5,t6) == CI(t5,t4)) ;
	            }
	          }
	        }
	      }
	    }
	  }
	}


	VERBOSE_DEBUG("pij = pij/iter;");

	pij = pij/iter;

	/**
	 * ## Compute Binder_f              ##
	 * ###################################
	 * Equivalent to:  for(g in 1:n_save){
	 *                     cij <- outer(c_out[[g]], c_out[[g]], "==");
	 *                     aux <- (pij - Const_Binder) * as.matrix(cij);
	 *                     aux <-  aux[upper.tri(aux)];
	 *                     Binder_f[g] <- sum(aux);
	 *                  }
	 *
	 */

	VERBOSE_DEBUG(" Compute Binder_f    ");

	arma::vec Binder_f (iter,arma::fill::zeros);

	arma::mat tmp = pij - Const_Binder;


	if ((N >= 1) && (iter >= 1)) {
	  for (t1=0;t1<=floord(iter-1,TILE_SIZE);t1++) {
	    for (t2=0;t2<=floord(N-1,TILE_SIZE);t2++) {
	      for (t3=TILE_SIZE*t2;t3<=std::min(N-1,TILE_SIZE*t2+(TILE_SIZE-1));t3++) {
	        for (t4=t3;t4<=N-1;t4++) {
	          lbv=TILE_SIZE*t1;
	          ubv=std::min(iter-1,TILE_SIZE*t1+(TILE_SIZE-1));
	          for (t5=lbv;t5<=ubv;t5++) {
	        	if (CI(t5,t4) == CI(t5,t3))
	        		Binder_f(t5) += tmp(t3,t4)   ;
	          }
	        }
	      }
	    }
	  }
	}




	/**
	 * ## Pick the maximum              ##
	 * ###################################
	 * Equivalent to:  	Binder_ind <- which.max(Binder_f)
	 *                  Binder_out <- c_out[Binder_ind]
	 *
	 */

	VERBOSE_DEBUG(" Return the maximum  ");
	arma::uword Binder_ind =  arma::index_max(Binder_f) ;


	return Binder_ind;

}


