#----------------------------------------------------------------------------
# Simple simulation example for binary network for PMP and MF variational inference
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
# Pearson correlation coefficient between true probabilities and estimated by MF
# Computational Cycles for PMP
# Parameter estimation coefficient between true probabilities and estimated by PMP

rm(list = ls())

set.seed(1234)




source('mean_field_DN.R')
source('helper.R')
source('mix_DN.R')

library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)



#--------------------------------Data Generation---------------------

n = 10   # n = 100; 

T = 10   # T = 100;

tau = 0.01 # tau =0.01 for the small transition case; tau =0.3 for the large transition case

rho = 0

beta = 0 # 

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = sqrt(0.5)  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(1.5,0,-1.5,0), nrow=2)

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
  Y[[t]] = positions_to_edges(X[[t]],beta)
}



#### Priors 

mean_beta_prior = 0    #prior mean of beta_0 

sigma_beta_prior = sqrt(2)   #prior sd of beta_0, 

sigma_X1 = sigma_0;   # initial sd of X[[1]], 

trans_sd = tau;   # sd of transition distribution, 

gap = 1e-2

Mix_list = mix_DN(Y = Y, rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior, sigma_X1 =sigma_0, trans_sd = tau ,
                  gap = gap)

MF_list =mean_field_DN (Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior, sigma_X1 =sigma_0, trans_sd = tau,
                          gap = gap)










# pearson correlation coeffcient
auc_mean_mf = rep((n-1)*n*T/2,0)
auc_mean_mix = rep((n-1)*n*T/2,0)
auc_res = rep((n-1)*n*T/2,0)
r=1
for(t in 1:T){
for (i in 1:n){
  for (j in 1:n){
    if (j<i){
      auc_mean_mf[r] = 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[T]][i,])%*% MF_list$Mean_X[[T]][j,]))
      auc_mean_mix[r] = 1/(1+exp(-Mix_list$mean_beta-t( Mix_list$Mean_X[[T]][i,])%*% Mix_list$Mean_X[[T]][j,]))
      auc_res[r] = 1/(1+exp(-beta-t(X[[T]][i,])%*% X[[T]][j,]))
      r=r+1
    }
  }
}
}
auc_mf <- cor(auc_res, auc_mean_mf,method = "pearson")
auc_mix <- cor(auc_res, auc_mean_mix,method = "pearson")

cat('Cycles for MF:',MF_list$iter,'\n')
cat('Pearson correlation for MF:',auc_mf,'\n')
cat('Cycles for PMP:',Mix_list$iter,'\n')
cat('Pearson correlation for PMP:',auc_mix,'\n')


