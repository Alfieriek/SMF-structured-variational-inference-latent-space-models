
rm(list = ls())

set.seed(2023)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(file_location)

source('mean_field_DN_adaptive.R')
source('helper.R')
source('mix_DN_adaptive.R')
source('mix_DN_adaptive_node.R')
source('mcmc_DN_adaptive.R')

library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(Bessel)


auc_node = NULL

auc_non_node = NULL

for (tau in c(0.05,0.1,0.5,1)){


  
for (simulation_id in 1:25){

#--------------------------------Data Generation---------------------

n = 20   # n = 100; 

T = 50   # T = 100;

#tau = 0.5 # tau=0.01,0.02,0.05,0.1,0.2,0.3

rho = 0

beta = 0 # 

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = sqrt(0.1)  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(1,0,-1,0), nrow=2)

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


for(i in sample(1:n,size=0.1*n)){
  eps =  rmvnorm(n=d,mean=rep(0,T), sigma= tau_Sigma)
  for (t in 2:T) {
    X[[t]][i,] <- X[[t-1]][i,] +eps[,t]
    #+  sample(c(-1,1),size=1)*0.01
  }
}



Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges(X[[t]],beta)
}



#### Priors 

mean_beta_prior = 0    #prior mean of beta_0 

sigma_beta_prior = sqrt(10)   #prior sd of beta_0, 


gap = 1e-3

start_time_MF <- Sys.time()


MF_list =mix_DN_adaptive (Y = Y,rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                           gap = gap)

end_time_MF <- Sys.time()




start_time_Mix <- Sys.time()


Mix_list = mix_DN_adaptive_node(Y = Y, rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                           gap = gap)

end_time_Mix <- Sys.time()





# pearson correlation coefficient
auc_mean_mf = rep((n-1)*n/2,0)
auc_mean_mix = rep((n-1)*n/2,0)
auc_res = rep((n-1)*n/2,0)

auc_mf = NULL


auc_mix = NULL

for(t in 1:T){
  r=1
  for (i in 1:n){
    for (j in 1:n){
      if (j<i){
        auc_mean_mf[r] = 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,]))
        auc_mean_mix[r] = 1/(1+exp(-Mix_list$mean_beta-t( Mix_list$Mean_X[[t]][i,])%*% Mix_list$Mean_X[[t]][j,]))
        auc_res[r] = 1/(1+exp(-beta-t(X[[t]][i,])%*% X[[t]][j,]))
        r=r+1
      }
    }
  }
  auc_mf <- c(auc_mf,cor(auc_res, auc_mean_mf,method = "pearson"))

  
  auc_mix <- c(auc_mix,cor(auc_res, auc_mean_mix,method = "pearson"))
}


cat('Cycles for non-nodewise SMF:',MF_list$iter,'\n')
cat('Running time for non-nodewise SMF:',end_time_MF-start_time_MF,'\n')
cat('Pearson correlation for non-nodewise SMF:',mean(auc_mf),'\n')


cat('Cycles for nodewise SMF:',Mix_list$iter,'\n')
cat('Running time for nodewise SMF:',end_time_Mix-start_time_Mix,'\n')
cat('Pearson correlation for nodewise SMF:',mean(auc_mix),'\n')


auc_node = c(auc_node,mean(auc_mix))
auc_non_node = c(auc_non_node,mean(auc_mf))
}
}
#boxplot(auc_node,auc_non_node)

tau_list = as.character(rep(c(rep(0.05,25),rep(0.1,25),rep(0.5,25),rep(1,25)),2))

 df = data.frame(list(tau = factor(tau_list),
                      method=c(rep('node-adaptive',100),rep('common-variance',100)),PCC = c(auc_node,auc_non_node)  ))
 

 library(ggplot2)
 ggplot(df, aes(x = tau, y = PCC, fill = method)) +
   geom_boxplot(outlier.shape = NA) +
   labs(x = "tau", y = "PCC scores", fill = "method",
        title = "Node-adaptive transition variance vs common transition variance") +
   theme_bw()+
   theme(plot.title = element_text(hjust = 0.5))
