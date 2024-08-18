#----------------------------------------------------------------------------
# Gaussian network via MCMC
# Input Variables:

# trans_sd: transition sd
# Y: network values, list of length T, Y[[t]]: n*n network values at time point t
# mean_beta_prior,sigma_beta_prior: mean and sd for prior of beta
# rho: ar1 coefficient, fixed to 1 in the paper
# sigma_X1: prior sd for the initial state
# sigma: sd for the likelihood
# gap: gap between errors for convergence
# max_iter: maximal cycles in computation
# alpha: fractional power of the likelihood, fixed to be 0.95 in this paper

#Output variables:
# err: dynamic of the training errors
# Mean_X: variational mean of each node, list of length T, Mean_X[[t]]: n*d matrix at time point t
# Sigma_X: variational covariance matrix of each node 
# mean_beta: variational mean of intercept beta
# sigma_beta: variational sd of intercept beta 
# iter: k-1


library(mvtnorm)


MCMC_GaussianDN = function(Y,trans_sd ,mean_beta_prior, sigma_beta_prior,rho=1, sigma_X1 ,sigma,  gap = 0.001, max_iter=10, alpha=0.95){
  
  T = length(Y)
  
  n = length(Y[[1]][1,])
  
  #### target parameters
  
  beta = 0 
  
  X = vector("list", T); 
  
  for (t in 1:T){
    X[[t]] = matrix(rep(0,n*d),nrow = n)
    for (i in 1:n){
      X[[t]][i,] =  rmvnorm(n=1,mean=rep(0,d), sigma=diag(d)) # randomly initialization for mu_{it}
    }
  }
  
  
  #--------------------------------Algorithm---------------------
  K= 500
  
  err =rep(0,K)
  
  err[1] =1e10
  
  k = 2
  
  ind = 0
  
  
  
  V_beta_cumulative = 0
  
  M_beta_cumulative = 0
  
  V_X_cumulative =  replicate(n=T, expr=list())
  
  M_X_cumulative = vector("list", T)
  
  
  for (t in 1:T){
    M_X_cumulative[[t]] = matrix(rep(0,n*d),nrow = n)
    V_X_cumulative[[t]] = vector("list", n)
    for (i in 1:n){
      V_X_cumulative[[t]][[i]] =  0*diag(d)        
    }
  }
  
  mean_beta_new = 0;    
  
  sigma_beta_new = 0  
  
  Mean_X_new = vector("list", T)
  
  Sigma_X_new = replicate(n=T, expr=list())
  
  
  
  for (t in 1:T){
    Mean_X_new[[t]] = matrix(rep(0,n*d),nrow = n)
    Sigma_X_new[[t]] = vector("list", n)
    for (i in 1:n){
      Sigma_X_new[[t]][[i]] =  0*diag(d)       
    }
  }
  
  
  while(k<K && ind ==0){
    
    
    # Calculate auxiliary values
    
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:n){
          M_i = X[[t]][i,]
          M_j = X[[t]][j,]
          if (j!= i){ 
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]- t(M_i) %*% M_j)*alpha
          }
        }
      }
    }
    
    
    
    ### update of sigma_beta
    
    sigma_beta_new = (sigma_beta_prior^(-2)+ T*n*(n-1)*sigma^(-2) )^(-1/2)
    
    
    ### update of mean_beta
    
    mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + sigma^(-2)*M_beta_cumulative) )
    
    
    beta = rnorm(n=1, mean= mean_beta_new,sd = sigma_beta_new)
    
    ### update of Sigma_X, Mean_X
    
    
    
    for (i in 1:n){
      
      for (t in 1:T){
        
        for (j in 1:n){
          
            M_i = X[[t]][i,]
            M_j = X[[t]][j,]

          if (j!= i){
            V_X_cumulative[[t]][[i]] = V_X_cumulative[[t]][[i]] +  sigma^(-2)*(M_j %*% t(M_j) ) *alpha
            M_X_cumulative[[t]][i,] = M_X_cumulative[[t]][i,] +  sigma^(-2)*(Y[[t]][i,j]-beta) * M_j *alpha
          }
        }
      
      
      
      if(t == 1)
      { Sigma_X_new[[t]][[i]] = ginv(sigma_X1^(-2)*diag(d)+trans_sd^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
      
      Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
      
      Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
        (trans_sd ^(-2)*X[[t+1]][i,]+M_X_cumulative[[t]][i,])
      }
      
      
      if(1<t  && t<T)
      { Sigma_X_new[[t]][[i]] = ginv(2*trans_sd ^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
      
      Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
      
      Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
        (trans_sd^(-2)*X[[t-1]][i,]+trans_sd ^(-2)*X[[t+1]][i,]+M_X_cumulative[[t]][i,])
      }
      if(t==T)
      { Sigma_X_new[[t]][[i]] = ginv(trans_sd ^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
      
      Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
      
      Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
        (trans_sd ^(-2)*X[[t-1]][i,]+M_X_cumulative[[t]][i,])
      }
        X[[t]][i,] = rmvnorm(n=1, mean = Mean_X_new[[t]][i,],sigma = Sigma_X_new[[t]][[i]])
      
      }  
      
    }
      
    pred_mean = rep(T*(n-1)*n,0)
    res = rep(T*(n-1)*n,0)
    r=1
    for (t in 1:T){
      for (i in 1:n){
        for (j in 1:n){
          if (j!=i){
            pred_mean[r] = mean_beta_new+t(Mean_X_new[[t]][i,])%*% Mean_X_new[[t]][j,]
            res[r] = Y[[t]][i,j]
            r=r+1
          }
        }
      }
    }
    
    
    err[k]=sqrt(sum((res-pred_mean)^2/T/n/(n-1)))
    
    if( k>50){
      ind = 1
    }    else{
      #err[k]-err[k-1] >-0.001
      
      # Mean_X = Mean_X_new
      # 
      # Sigma_X = Sigma_X_new
      # 
      # mean_beta = as.numeric( mean_beta_new)
      # 
      # sigma_beta = sigma_beta_new
      # 
      
      cat(err[k],'\n')
      
      
      cat('iteration:',k-1,'\n')
      
      
      k = k+1
    }
    
  }
  
  return(list(err=err[2:(k-1)], Mean_X= Mean_X_new, Sigma_X= Sigma_X_new,  mean_beta  = mean_beta_new, sigma_beta =sigma_beta_new , iter = k-1 ))
}
