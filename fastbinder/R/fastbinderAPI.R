
#################################################################################
##### Binder Functions
#################################################################################

#'  Run the naive binder algorithm (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder_naive=function(CI,  Const_Binder = 0.5){
	return (IAM_binder_naive(CI, Const_Binder))
}

#'  Run the ptimized binder algorithm (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder_opt=function(CI,  Const_Binder = 0.5){
	return (IAM_binder_opt(CI, Const_Binder))
}

#'  Run the parallel binder algorithm (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder_parallel=function(CI,  Const_Binder = 0.5){
	return (IAM_binder_parallel(CI, Const_Binder))
}

#'  Run the parallel binder algorithm (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder = AM_binder_parallel

#'  Run the binder algorithm using Andrea implementation (TBD)
#'  
#'  TBD
#'  
#'@param C fit CI
#'@param Const_Binder Const_Binder
#'  
#'@export
AM_binder_andrea=function (CI,  weight = 0.5) {
	
	#Equal costs
	Const_Binder <- weight
	
	#c_out contains matrix of allocation labels
	c_out <- CI
	
	
	#c_contains matrix of allocation labels
	
	
	n_save <- dim(c_out)[1]
	N      <- dim(c_out)[2]
	
	
	
	#Compute similarity matrix pij
	
	pij <- matrix(0,N,N)
	
	for(g in 1:n_save){
		
		pij <- pij + outer(c_out[g,], c_out[g,], "==")
		
	}
	
	pij <- pij/n_save
	
	
	
	Binder_f <- rep(0,n_save)
	
	for(g in 1:n_save){
		
		cij <- outer(c_out[g,], c_out[g,], "==")
		
		aux <- (pij - Const_Binder) * as.matrix(cij)
		
		aux <-  aux[upper.tri(aux)]
		
		Binder_f[g] <- sum(aux)
		
	}
	
	Binder_ind <- which.max(Binder_f)
	
	
	
	return (Binder_ind)
}
