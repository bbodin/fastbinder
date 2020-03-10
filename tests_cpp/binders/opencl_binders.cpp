/*
 * opencl_binders.cpp
 *
 *  Created on: Mar 9, 2020
 *      Author: toky
 */



#include "binders.h"


#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif



/**
 * Naive implementation as a baseline
 */
arma::uword  opencl_binder (arma::mat CI , double Const_Binder ) {

  


	  // Get platform and device information
	  cl_platform_id * platforms = NULL;
	  cl_uint     num_platforms;
	  //Set up the Platform
	  cl_int clStatus = clGetPlatformIDs(0, NULL, &num_platforms);
	  platforms = (cl_platform_id *)
	  malloc(sizeof(cl_platform_id)*num_platforms);
	  clStatus = clGetPlatformIDs(num_platforms, platforms, NULL);

	  //Get the devices list and choose the device you want to run on
	  cl_device_id     *device_list = NULL;
	  cl_uint           num_devices;

	  clStatus = clGetDeviceIDs( platforms[0], CL_DEVICE_TYPE_GPU, 0,NULL, &num_devices);
	  device_list = (cl_device_id *)
	  malloc(sizeof(cl_device_id)*num_devices);
	  clStatus = clGetDeviceIDs( platforms[0],CL_DEVICE_TYPE_GPU, num_devices, device_list, NULL);

	  // Create one OpenCL context for each device in the platform
	  cl_context context;
	  context = clCreateContext( NULL, num_devices, device_list, NULL, NULL, &clStatus);

	  // Create a command queue
	  cl_command_queue command_queue = clCreateCommandQueue(context, device_list[0], 0, &clStatus);





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


