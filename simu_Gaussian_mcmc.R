#----------------------------------------------------------------------------
# Simple simulation example for Gaussian network for SMF and MF variational inference
# Key Variables:
# n: number of nodes
# T: number of time points
# tau: transition sd
# X: true latent positions, list of length T, X[[t]]: n*d latent matrix at time point t
# Y: network values, list of length T, Y[[t]]: n*n network values at time point t
# MF_list: model trained by MF 
# Mix_list: model trained by SMF

#Output Comparisons:
# Computational Cycles for MF 
# Parameter estimation error for MF
# Computational Cycles for SMF
# Parameter estimation error for SMF


rm(list = ls())
set.seed(2023)
library(MASS)
library(mvtnorm)
library(igraph)
library(Bessel)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(file_location)
source('helper.R')
source('mean_field_Gaussian_adaptive.R')
source('mix_Gaussian_adaptive.R')
source('mcmc_Gaussian_adaptive.R')
#source('mix_Gaussian_adaptive_node.R')
# 

#--------------------------------Data Generation---------------------

n = 20   # n = 100;  

T = 100   # T = 100;

tau = 0.1

rho = 0.5 

beta = -1 # 


d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = tau  # sd for initial state

sigma = 0.1 #sd for likelihood of y

 initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

 initial_mean <- matrix(c(0,0,0,0), nrow=2)

initial_Sigma <- sigma_0^2*diag(d)

X <- vector("list", T)

for ( i in 1:n){
  X[[1]] = rbind(X[[1]], rmvnorm(n=1,mean=initial_mean[,initial_components[i]], sigma=initial_Sigma))
}

for(t in 1:(T-1)){X[[t+1]] = X[[1]]}


sig_eps = diag(rep(1,T))

for ( t1 in 1:(T)){
  for (t2 in 1:(T)){
    if (abs(t1-t2)!=0){sig_eps[t1,t2] = rho }
  }
}

tau_Sigma = tau^2*sig_eps


for(i in 1:n){
  eps =  rmvnorm(n=d,mean=rep(0,T), sigma= tau_Sigma)
  for (t in 2:T) {
    X[[t]][i,] <- X[[t-1]][i,] +  eps[,t-1]
  }
}


Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_Gaussian(X[[t]],beta,sigma)
}



#--------------------------------Model Fitting---------------------


mean_beta_prior = 0    #prior mean of beta, 

sigma_beta_prior = 10   #prior sd of beta, 

sigma_X1 = sigma_0;   # initial sd of X[[1]], 

trans_sd = tau;   # sd of transition distribution, 

gap =1e-4

min_iter = 1


start_time_MF <- Sys.time()

MF_list_adaptive = mcmc_GaussianDN_adaptive(Y=Y ,mean_beta_prior=mean_beta_prior, d=d,
                                            max_iter=100, sigma_beta_prior=sigma_beta_prior,sigma= sigma,gap= 1e-4)
end_time_MF <- Sys.time()

# New code for max_iter=50
start_time_MF_50 <- Sys.time()

MF_list_adaptive_50 = mcmc_GaussianDN_adaptive(Y=Y ,mean_beta_prior=mean_beta_prior, d=d,
                                               max_iter=50, sigma_beta_prior=sigma_beta_prior,sigma= sigma,  gap= 1e-4)
end_time_MF_50 <- Sys.time()

# New code for max_iter=200
start_time_MF_200 <- Sys.time()

MF_list_adaptive_200 = mcmc_GaussianDN_adaptive(Y=Y ,mean_beta_prior=mean_beta_prior, d=d,
                                                max_iter=200, sigma_beta_prior=sigma_beta_prior,sigma= sigma,  gap= 1e-4)
end_time_MF_200 <- Sys.time()


start_time_Mix <- Sys.time()

Mix_list_adaptive = mix_Gaussian_adaptive(Y=Y ,mean_beta_prior=mean_beta_prior, sigma_beta_prior=sigma_beta_prior,
                                          global_prior='Cauthy',sigma = sigma, gap = 1e-4)

end_time_Mix <- Sys.time()



pred_mean_MF_adaptive = matrix(rep(0,T*n*(n-1)),nrow = T)
pred_mean_MF_adaptive_50 = matrix(rep(0,T*n*(n-1)),nrow = T)
pred_mean_MF_adaptive_200 = matrix(rep(0,T*n*(n-1)),nrow = T)

pred_mean_Mix_adaptive = matrix(rep(0,T*n*(n-1)),nrow = T)
res =  matrix(rep(0,T*n*(n-1)),nrow = T)


MF_adaptive_sm = rep(0,T)
Mix_adaptive_sm = rep(0,T)
MF_adaptive_sm_50 = rep(0,T)
MF_adaptive_sm_200 = rep(0,T)


for (t in 1:T){
  r=1
  for (i in 1:n){
    for (j in 1:n){
      if (j!=i){
        # Compute predicted means for new MF cases
        pred_mean_MF_adaptive_50[t,r] = t(MF_list_adaptive_50$Mean_X[[t]][i,])%*% MF_list_adaptive_50$Mean_X[[t]][j,]+MF_list_adaptive_50$mean_beta
        pred_mean_MF_adaptive_200[t,r] = t(MF_list_adaptive_200$Mean_X[[t]][i,])%*% MF_list_adaptive_200$Mean_X[[t]][j,]+MF_list_adaptive_200$mean_beta
        pred_mean_MF_adaptive[t,r] = t(MF_list_adaptive$Mean_X[[t]][i,])%*% MF_list_adaptive$Mean_X[[t]][j,]+MF_list_adaptive$mean_beta
        pred_mean_Mix_adaptive[t,r] = t(Mix_list_adaptive$Mean_X[[t]][i,])%*% Mix_list_adaptive$Mean_X[[t]][j,]+Mix_list_adaptive$mean_beta
        res[t,r] = t(X[[t]][i,]) %*% X[[t]][j,]+beta
        r=r+1
      }
    }
  }
  MF_adaptive_sm_50[t] = sum((pred_mean_MF_adaptive_50[t,]-res[t,])^2/n/T/(n-1))
  MF_adaptive_sm_200[t] = sum((pred_mean_MF_adaptive_200[t,]-res[t,])^2/n/T/(n-1))
  MF_adaptive_sm[t] = sum((pred_mean_MF_adaptive[t,]-res[t,])^2/n/T/(n-1))
  Mix_adaptive_sm[t] = sum((pred_mean_Mix_adaptive[t,]-res[t,])^2/n/T/(n-1))
}





cat('Cycles for adaptive SMF:',Mix_list_adaptive$iter,'\n')
cat('Running time for SMF:',end_time_Mix-start_time_Mix,'\n')
cat('Parameter estimation error for adaptive SMF:',sqrt(sum(Mix_adaptive_sm)),'\n')


cat('Cycles for adaptive MF with 50 iterations:',MF_list_adaptive_50$iter,'\n')
cat('Running time for adaptive MF with 50 iterations:',end_time_MF_50-start_time_MF_50,'\n')
cat('Parameter estimation error for adaptive MF with 50 iterations:',sqrt(sum(MF_adaptive_sm_50)),'\n')


cat('Cycles for adaptive MF with 100 iterations:',MF_list_adaptive$iter,'\n')
cat('Running time for adaptive MF:',end_time_MF-start_time_MF,'\n')
cat('Parameter estimation error for adaptive MF:',sqrt(sum(MF_adaptive_sm)),'\n')


cat('Cycles for adaptive MF with 200 iterations:',MF_list_adaptive_200$iter,'\n')
cat('Running time for adaptive MF with 200 iterations:',as.numeric(end_time_MF_200-start_time_MF_200,units = "secs"),'\n')
cat('Parameter estimation error for adaptive MF with 200 iterations:',sqrt(sum(MF_adaptive_sm_200)),'\n')


