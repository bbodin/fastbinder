quit() ## Skip this test - too long
##############################################
### Load the AntMan package
##############################################

rm(list=ls())
.libPaths(c("/home/toky/yalenus/research/mixture/AntMan/AntMAN.Rinstall", .libPaths()))
library("AntMAN")

##############################################
### BUILD THE MULTIVARIATE NORMAL DATA
##############################################

set.seed(123)


MU <- matrix(nrow=3,ncol=2)

MU[1,] <- c(0,0)
MU[2,] <- c(-3,-3)
MU[3,] <- c(4,4)


sig1 <- c(1,1)
rho1 <- 0
Sig1 <- matrix(c(sig1[1]^2,rho1*sig1[1]*sig1[2], rho1*sig1[1]*sig1[2],sig1[2]^2),byrow=TRUE,nrow=2) 

sig2 <- c(1,1)
rho2 <- -0.7
Sig2 <- matrix(c(sig2[1]^2,rho2*sig2[1]*sig2[2], rho2*sig2[1]*sig2[2],sig2[2]^2),byrow=TRUE,nrow=2) 

sig3 <- c(1,1)
rho3 <- -0.3
Sig3 <- matrix(c(sig3[1]^2,rho3*sig3[1]*sig3[2], rho3*sig3[1]*sig3[2],sig3[2]^2),byrow=TRUE,nrow=2) 


SIG <- array(0,dim=c(3,2,2))
SIG[1,,] <- Sig1
SIG[2,,] <- Sig2
SIG[3,,] <- Sig3


N = 1000

demo_multivariate_normal <-AM_sample_multinorm(n = N ,d = 2,c(0.3,0.3,0.4),MU,SIG)
y_mvn  <- demo_multivariate_normal$y
ci_mvn <- demo_multivariate_normal$ci

hist(y_mvn,freq=FALSE,nclass=15,col=colors()[4])
plot(y_mvn,col=ci_mvn+1)


##############################################################################
### PREPARE THE GIBBS for multivariate Normal mixture with poisson gamma priors
##############################################################################


mixture_mvn_params = AM_multinorm_mix_hyperparams   (mu0=c(0,0),ka0=1,nu0=4,Lam0=diag(2))

mcmc_params        = AM_mcmc_parameters(niter=4000, burnin=1500, thin=1, verbose=1, output=c("CI"))
components_prior   = AM_mix_components_prior_pois (init=3,  a=1, b=1) 
weights_prior      = AM_mix_weights_prior_gamma(init=2, a=1, b=1)

fit <- AM_mcmc_fit(
  y = y_mvn, 
  mix_kernel_hyperparams = mixture_mvn_params,
  mix_components_prior =components_prior,
  mix_weight_prior = weights_prior,
  mcmc_parameters = mcmc_params)



nrows = length(fit$CI)
ncols = length(fit$CI[[1]])

data = array(as.numeric(unlist(fit$CI)), dim=c(nrows, ncols))

dim(data)

##############################################################################
### Run Andrea's
##############################################################################




start_time1 <- Sys.time()
AM_binder_andrea(data)
end_time1 <- Sys.time()

start_time2 <- Sys.time()
AM_binder_naive(data)
end_time2 <- Sys.time()

start_time3 <- Sys.time()
AM_binder_opt(data)
end_time3 <- Sys.time()


start_time4 <- Sys.time()
AM_binder_parallel(data)
end_time4 <- Sys.time()


end_time1 - start_time1
end_time2 - start_time2
end_time3 - start_time3
end_time4 - start_time4

