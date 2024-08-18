rm(list = ls())

set.seed(1234)


file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(file_location)

source('mean_field_DN_adaptive.R')
source('helper.R')
source('mix_DN_adaptive.R')

library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)



#--------------------------------Data Generation---------------------

n = 20   # n = 100; 

T = 100   # T = 100;

tau = 0.01 # tau =0.01 for the small transition case; tau =0.3 for the large transition case

rho = 0

beta = 0 # 

d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = sqrt(0.1)  # sd for initial state

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(0.5,0,-0.5,0), nrow=2)

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

gap = 1e-3

Mix_list = mix_DN_adaptive(Y = Y, rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior, 
                  gap = gap)

Mix_list_beta_zero =mix_DN_adaptive(Y = Y, rho=1, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior, 
                                    gap = gap, beta_zero = TRUE)









# pearson correlation coeffcient
auc_mean_mix_beta_zero = rep((n-1)*n*T/2,0)
auc_mean_mix = rep((n-1)*n*T/2,0)
auc_res = rep((n-1)*n*T/2,0)
r=1


for(t in 1:T){
for (i in 1:n){
  for (j in 1:n){
    if (j<i){
      auc_mean_mix_beta_zero[r] = -t( Mix_list_beta_zero$Mean_X[[T]][i,])%*% Mix_list_beta_zero$Mean_X[[T]][j,]
      auc_mean_mix[r] = -t( Mix_list$Mean_X[[T]][i,])%*% Mix_list$Mean_X[[T]][j,]
      auc_res[r] = -t(X[[T]][i,])%*% X[[T]][j,]
      r=r+1
    }
  }
}
}
auc_mix_beta_zero <- cor(auc_res, auc_mean_mix_beta_zero,method = "pearson")
auc_mix <- cor(auc_res, auc_mean_mix,method = "pearson")


par(mfrow=c(1,2))


xli = max(abs(auc_res))
yli = max(abs(c(auc_mean_mix)))
plot(auc_res,auc_mean_mix,main="Unknown beta",xlab="True",
     ylab="Est",xlim = c(-xli,xli),ylim=c(-yli,yli))
lines(seq(from=-xli,to=xli,by=xli/50),seq(from=-xli,to=xli,by=xli/50),col=2)

xli = max(abs(auc_res))
yli = max(abs(c(auc_mean_mix_beta_zero)))
plot(auc_res,auc_mean_mix,main="Known beta",xlab="True",
     ylab="Est",xlim = c(-xli,xli),ylim=c(-yli,yli))
lines(seq(from=-xli,to=xli,by=xli/50),seq(from=-xli,to=xli,by=xli/50),col=2)

cat("The estimated beta for the unknown beta is", Mix_list$mean_beta, "\n")
