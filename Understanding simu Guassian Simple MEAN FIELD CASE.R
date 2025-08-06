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
set.seed(2020)
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
set.seed(2020)
n = 4   # n = 100;  

T = 5   # T = 100;

tau = 1.2 # tau=0.01,0.02,0.05,0.1,0.2,0.3,0.5

rho = 0 

beta = 0 # 


d = 2   # dimension for the latent vectors, d = 2 for better visualization

sigma_0 = tau  # sd for initial state

sigma = 0.1 #sd for likelihood of y

initial_components = sample(1:2,prob=c(0.5,0.5),size=n, replace=TRUE)

initial_mean <- matrix(c(0,0,0,0), nrow=2)

initial_Sigma <- sigma_0^2*diag(d)

X <- vector("list", T)

# This is essentially sampling 10 points from the circle with radius 0.0001.
for ( i in 1:n){
  X[[1]] = rbind(X[[1]], rmvnorm(n=1,mean=initial_mean[,initial_components[i]], sigma=initial_Sigma))
}

# Plotting First Latent Position

# plot(X[[1]][,1],X[[1]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1), main = paste("Time Step",1),pch = 16,cex = 1)
# text(X[[1]][,1], X[[1]][,2], c(1:n),cex = 0.5, pos = 2, col = "blue")
# abline(h = 0, v = 0, col = "black", lty = 2)

# Add contour lines
# z <- outer(seq(-0.1, 0.1, length.out = 100), seq(-0.1, 0.1, length.out = 100), function(x, y) {
#   dmvnorm(cbind(x, y), mean = initial_mean[,initial_components[i]], sigma = initial_Sigma)
# })
# contour(seq(-0.1, 0.1, length.out = 100), seq(-0.1, 0.1, length.out = 100), z , add = T, drawlabels = F, col = "red")



# Just populating every entries in the list if 1:T elements with X[[1]] for all. 
for(t in 1:(T-1)){X[[t+1]] = X[[1]]}


sig_eps = diag(rep(1,T))

for ( t1 in 1:(T)){
  for (t2 in 1:(T)){
    if (abs(t1-t2)!=0){sig_eps[t1,t2] = rho }
  }
}

tau_Sigma = tau^2*sig_eps

# Updating the ith row of the X[[t]] for t = 2:T in the loop i = 1:n.
for(i in 1:n){
  eps =  rmvnorm(n=d,mean=rep(0,T), sigma= tau_Sigma)
  for (t in 2:T) {
    X[[t]][i,] <- X[[t-1]][i,] +  eps[,t-1]
  }
}

# Plotting Latent Positions
par(mfrow = c(4,2))
for (i in 1:T) {
  plot(X[[i]][,1],X[[i]][,2], xlim = c(min(unlist(X)),max(unlist(X))), ylim = c(min(unlist(X)),max(unlist(X))), main = paste("Time Step",i),pch = 16,cex = 1.5)
  text(X[[i]][,1], X[[i]][,2], c(1:n),cex = 1, pos = 2, col = "blue")
  abline(h = 0, v = 0, col = "black", lty = 2)
}



Y = vector("list", T)

for (t in 1:T){
  Y[[t]] = positions_to_edges_Gaussian(X[[t]],beta,sigma)
}


# positions_to_edges_Gaussian simple example
# A <- X[[1]]
#   n = length(A[,1])
#   edge_mat = matrix(rep(0,n*n),nrow=n)
# 
#   for (i in 1:n){
#     for(j in 1:n){
#       if(j > i)
#       {edge_mat[j,i] = rnorm(n=1,mean=beta+t(A[i,])%*%A[j,],sd=sigma)}
#     }
#   }
#   edge_mat = edge_mat+ t(edge_mat)
# 



#--------------------------------Model Fitting--------------------- MEAN FIELD


mean_beta_prior = 0    #prior mean of beta, 

sigma_beta_prior = 10   #prior sd of beta, 

sigma_X1 = sigma_0;   # initial sd of X[[1]], 

trans_sd = tau;   # sd of transition distribution, 

gap =1e-3

min_iter = 1


# Running the Mean Field Variational inference for dynamic latent space model

MF_SMF_EG <- mean_field_GaussainDN_adaptive(Y=Y,
                                            mean_beta_prior= mean_beta_prior,
                                            sigma_beta_prior = sigma_beta_prior,
                                            rho=1,
                                            sigma = sigma ,
                                            gap =0.01,
                                            max_iter=100,
                                            alpha=0.95)


MF_SMF_EG_MEANs <- MF_SMF_EG$Mean_X

# Plotting Latent Positions
par(mfrow = c(4,2))
for (i in 1:T) { 
  plot(X[[i]][,1],X[[i]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1), main = paste("Time Step",i),pch = 16,cex = 1.5, col = "blue")
  text(X[[i]][,1], X[[i]][,2], c(1:n),cex = 1, pos = 2, col = "blue")
  abline(h = 0, v = 0, col = "black", lty = 2)
  points(MF_SMF_EG_MEANs[[i]][,1],MF_SMF_EG_MEANs[[i]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1),pch = 18,cex = 1.5, col = "red")
  text(MF_SMF_EG_MEANs[[i]][,1], MF_SMF_EG_MEANs[[i]][,2], c(1:n),cex = 1, pos = 2, col = "red")
  legend("topleft", legend=c("TRUE","MF_VI"), col=c("blue", "red"), pch=c(16,18), cex=0.55)
}

# Running the Structured Mean Field Variational inference for dynamic latent space model

SMF_SMF_EG <- mix_Gaussian_adaptive(Y=Y ,
                      mean_beta_prior=mean_beta_prior,
                      sigma_beta_prior=sigma_beta_prior,
                      sigma = sigma,
                      X_true=X,
                      gap = 1e-6)


SMF_SMF_EG_MEANs <- SMF_SMF_EG$Mean_X

# Plotting Latent Positions
# par(mfrow = c(2,2))
# for (i in 1:T) { 
#   plot(X[[i]][,1],X[[i]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1), main = paste("Time Step",i),pch = 16,cex = 1.5, col = "blue")
#   text(X[[i]][,1], X[[i]][,2], c(1:n),cex = 1, pos = 2, col = "blue")
#   abline(h = 0, v = 0, col = "black", lty = 2)
#   points(MF_SMF_EG_MEANs[[i]][,1],MF_SMF_EG_MEANs[[i]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1),pch = 18,cex = 1.5, col = "red")
#   text(MF_SMF_EG_MEANs[[i]][,1], MF_SMF_EG_MEANs[[i]][,2], c(1:n),cex = 1, pos = 2, col = "red")
#   points(SMF_SMF_EG_MEANs[[i]][,1],SMF_SMF_EG_MEANs[[i]][,2], xlim = c(-0.1,0.1), ylim = c(-0.1,0.1),pch = 20,cex = 1.5, col = "orange")
#   text(SMF_SMF_EG_MEANs[[i]][,1], SMF_SMF_EG_MEANs[[i]][,2], c(1:n),cex = 1, pos = 2, col = "orange")
#   legend("topleft", legend=c("TRUE","MF_VI"), col=c("blue", "red","orange"), pch=c(16,18,20), cex=0.55)
# }




# Plotting Latent Positions
par(mfrow = c(2,1))
for (i in 1:T) { 
  plot(X[[i]][,1],X[[i]][,2], xlim = c(min(X[[i]][,1],MF_SMF_EG_MEANs[[i]][,1],SMF_SMF_EG_MEANs[[i]][,1]),
                                       max(X[[i]][,1],MF_SMF_EG_MEANs[[i]][,1],SMF_SMF_EG_MEANs[[i]][,1])),
                                       ylim = c(min(X[[i]][,2],MF_SMF_EG_MEANs[[i]][,2],SMF_SMF_EG_MEANs[[i]][,2]),
                                                max(X[[i]][,2],MF_SMF_EG_MEANs[[i]][,2],SMF_SMF_EG_MEANs[[i]][,2])),
                                                main = paste("Time Step",i),pch = 16,cex = 1.5, col = "blue")
  text(X[[i]][,1], X[[i]][,2], c(1:n),cex = 1, pos = 2, col = "blue")
  abline(h = 0, v = 0, col = "black", lty = 2)
  points(MF_SMF_EG_MEANs[[i]][,1],MF_SMF_EG_MEANs[[i]][,2],pch = 18,cex = 1.5, col = "red")
  text(MF_SMF_EG_MEANs[[i]][,1], MF_SMF_EG_MEANs[[i]][,2], c(1:n),cex = 1, pos = 2, col = "red")
  points(SMF_SMF_EG_MEANs[[i]][,1],SMF_SMF_EG_MEANs[[i]][,2],pch = 10,cex = 1.5, col = "orange")
  text(SMF_SMF_EG_MEANs[[i]][,1], SMF_SMF_EG_MEANs[[i]][,2], c(1:n),cex = 1, pos = 2, col = "orange")
  legend("topright", legend=c("TRUE","MF_VI","SMF_VI"), col=c("blue", "red","orange"), pch=c(16,18,10), cex=0.55)
}






