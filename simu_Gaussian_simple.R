#----------------------------------------------------------------------------
# Simple simulation example for Gaussian network for PMP and MF variational inference
# Key Variables:
# n: number of nodes
# T: number of time points
# tau: transition sd
# X: true latent positions, list of length T, X[[t]]: n*d latent matrix at time point t
# Y: network values, list of length T, Y[[t]]: n*n network values at time point t
# MF_list: model trained by MF 
# Mix_list: model trained by PMP

#Output Comparisons:
# Computational Cycles for MF 
# Parameter estimation error for MF
# Computational Cycles for PMP
# Parameter estimation error for PMP


rm(list = ls())
set.seed(2020)
library(MASS)
library(mvtnorm)
library(igraph)



source('mean_field_GaussianDN.R')
source('mix_GaussianDN.R')
source('helper.R')
# 

#--------------------------------Data Generation---------------------

n = 10   # n = 100;  

T = 10   # T = 100;

tau = 0.01   # tau=0.01/ tau=0.3

rho = 0 

beta = 0 # 


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

sigma_beta_prior = sqrt(2)   #prior sd of beta, 

sigma_X1 = sigma_0;   # initial sd of X[[1]], 

trans_sd = tau;   # sd of transition distribution, 

gap =1e-3

min_iter = 1

MF_list = mean_field_GaussainDN (Y=Y,trans_sd=trans_sd ,mean_beta_prior=mean_beta_prior, sigma_beta_prior=sigma_beta_prior,sigma= sigma,rho=1, sigma_X1 =sigma_X1,  gap= gap)

Mix_list = mix_GaussianDN (Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior, sigma_X1 =sigma_0, trans_sd = trans_sd,sigma =sigma,
                           gap = gap)

#--------------------------------Computation Efficiency and Estimation Error Comparison---------------------

pred_mean_MF = matrix(rep(0,T*n*(n-1)),nrow = T)
pred_mean_Mix = matrix(rep(0,T*n*(n-1)),nrow = T)
res =  matrix(rep(0,T*n*(n-1)),nrow = T)

MF_rmse_sm = rep(0,T)
Mix_rmse_sm = rep(0,T)

for (t in 1:T){
  r=1
for (i in 1:n){
  for (j in 1:n){
    if (j!=i){
      pred_mean_MF[t,r] = t(MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,]
      pred_mean_Mix[t,r] = t(Mix_list$Mean_X[[t]][i,])%*% Mix_list$Mean_X[[t]][j,]
      res[t,r] = beta+t(X[[t]][i,]) %*% X[[t]][j,]
      r=r+1
    }
  }
}
  MF_rmse_sm[t] = sum((pred_mean_MF[t,]-res[t,])^2/n/T/(n-1))
  Mix_rmse_sm[t] = sum((pred_mean_Mix[t,]-res[t,])^2/n/T/(n-1))
}

cat('Cycles for MF:',MF_list$iter,'\n')
cat('Parameter estimation error for MF:',sqrt(sum(MF_rmse_sm)),'\n')

cat('Cycles for PMP:',Mix_list$iter,'\n')
cat('Parameter estimation error for PMP:',sqrt(sum(Mix_rmse_sm)),'\n')



