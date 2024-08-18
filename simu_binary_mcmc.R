rm(list = ls())

set.seed(1234)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(file_location)

source('mean_field_DN_adaptive.R')
source('helper.R')
source('mix_DN_adaptive.R')
#source('mix_DN_adaptive_node.R')
source('mcmc_DN_adaptive.R')

library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)
library(BayesLogit)

n= 10

T = 100

prob = 0.95

beta = 0 

d = 2   

sigma = 0.3


X <- vector("list", n)

trend_generate = function(T,d,sigma_0,prob){
  X = matrix(rep(0,T*d),nrow = T)
  initial_components = sample(1:2,prob=c(0.5,0.5),size=1, replace=TRUE)
  initial_mean <- matrix(c(1,0,-1,0), nrow=2)
  X[1,] = rnorm(n=1, mean = initial_mean[,initial_components], sd = rep(sigma_0,d) )
  rho = prob
  for (t in 2:T) {
    increment = rep(sample(c(0,0.5),prob = c(rho, 1-rho),size=1),d)
    X[t,] <- X[t-1,]+increment*sign(runif(d,-1,1))
  }
  return(X)
}


for ( i in 1:n){
  X[[i]] = trend_generate(T,d,sigma,prob)
  
}



positions_to_edges_binary = function(X,beta,sigma,t){
  n = length(X)
  edge_mat = matrix(rep(0,n*n),nrow = n)
  for (i in 1:n){
    for(j in 1:n){
      prob = 1/(1+exp(-beta-t(X[[i]][t,])%*%X[[j]][t,]))
      edge_mat[i,j] = rbinom(n=1, size=1, prob=prob)
    }
  }
  return(edge_mat)
}


Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_binary(X,beta,sigma,t)
}


#### Priors 

mean_beta_prior = 0    #prior mean of beta_0 

sigma_beta_prior = sqrt(10)   #prior sd of beta_0, 


gap = 1e-4

start_time_MF <- Sys.time()


MF_list =mcmc_DN_adaptive (Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                          gap = gap,max_iter=100)

end_time_MF <- Sys.time()


# New code for max_iter=50
start_time_MF_50 <- Sys.time()

MF_list_50 = mcmc_DN_adaptive(Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                                               gap = gap,max_iter=50)
end_time_MF_50 <- Sys.time()

# New code for max_iter=200
start_time_MF_200 <- Sys.time()

MF_list_200 = mcmc_DN_adaptive(Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                                                gap = gap,max_iter=200)
end_time_MF_200 <- Sys.time()



start_time_Mix <- Sys.time()


Mix_list = mix_DN_adaptive(Y = Y, rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                           gap = gap)

end_time_Mix <- Sys.time()





# pearson correlation coeffcient
auc_mean_mf = rep((n-1)*n/2,0)
auc_mean_mf_50 = rep((n-1)*n/2,0)
auc_mean_mf_200 = rep((n-1)*n/2,0)
auc_mean_mix = rep((n-1)*n/2,0)
auc_res = rep((n-1)*n/2,0)

auc_mf = NULL
auc_mf_50 = NULL
auc_mf_200 = NULL

auc_mix = NULL

for(t in 1:T){
  r=1
for (i in 1:n){
  for (j in 1:n){
    if (j<i){
      auc_mean_mf[r] = 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,]))
      auc_mean_mf_50[r] = 1/(1+exp(-MF_list_50$mean_beta-t( MF_list_50$Mean_X[[t]][i,])%*% MF_list_50$Mean_X[[t]][j,]))
      auc_mean_mf_200[r] = 1/(1+exp(-MF_list_200$mean_beta-t( MF_list_200$Mean_X[[t]][i,])%*% MF_list_200$Mean_X[[t]][j,]))  
      auc_mean_mix[r] = 1/(1+exp(-Mix_list$mean_beta-t( Mix_list$Mean_X[[t]][i,])%*% Mix_list$Mean_X[[t]][j,]))
      auc_res[r] = 1/(1+exp(-beta-t(X[[t]][i,])%*% X[[t]][j,]))
      r=r+1
    }
  }
}
  auc_mf <- c(auc_mf,cor(auc_res, auc_mean_mf,method = "pearson"))
  auc_mf_50 <- c(auc_mf_50,cor(auc_res, auc_mean_mf_50,method = "pearson"))
  auc_mf_200 <- c(auc_mf_200,cor(auc_res, auc_mean_mf_200,method = "pearson"))
  
  auc_mix <- c(auc_mix,cor(auc_res, auc_mean_mix,method = "pearson"))
}

cat('Cycles for MF_50:',MF_list_50$iter,'\n')
cat('Running time for adaptive MF_50:',end_time_MF_50-start_time_MF_50,'\n')
cat('Pearson correlation for MF_50:',mean(auc_mf_50),'\n')

cat('Cycles for MF:',MF_list$iter,'\n')
cat('Running time for adaptive MF:',end_time_MF-start_time_MF,'\n')
cat('Pearson correlation for MF:',mean(auc_mf),'\n')

cat('Cycles for MF_200:',MF_list_200$iter,'\n')
cat('Running time for adaptive MF_200:',end_time_MF_200-start_time_MF_200,'\n')
cat('Pearson correlation for MF_200:',mean(auc_mf_200),'\n')

cat('Cycles for SMF:',Mix_list$iter,'\n')
cat('Running time for SMF:',end_time_Mix-start_time_Mix,'\n')

cat('Pearson correlation for SMF:',mean(auc_mix),'\n')

