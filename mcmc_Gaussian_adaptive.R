
#----------------------------------------------------------------------------
# Gaussian network via MF variational inference
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

library(GIGrvg)
library(invgamma)

mcmc_GaussianDN_adaptive = function(Y,mean_beta_prior, d=2,sigma_beta_prior,rho=1, X_true=NULL, sigma , gap =0.01,max_iter=100, alpha=0.95){
  
  
  
  T = length(Y)
  
  n = length(Y[[1]][1,])
  
  
  #--------------------------------Model Initialization---------------------
  
  #### target parameters
  
  
  
  mean_beta = 0;    #mean of beta, 
  
  sigma_beta = sqrt(2);  #sd of beta, 
  
  Mean_X = vector("list", T); #Mean of X, Mean_X[[t]][i,] 
  
  X_out = vector("list", T)
  
  for (t in 1:T){
    Mean_X[[t]] = matrix(rep(0,n*d),nrow = n)
    X_out[[t]] = matrix(rep(0,n*d),nrow = n)
    for (i in 1:n){
      Mean_X[[t]][i,] =  rmvnorm(n=1,mean=rep(0,d), sigma=diag(d)) # randomly initialization for mu_{it}
    }
  }
  
 beta_out = 0
  
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
  
  #--------------------------------Algorithm---------------------
  K= 1000
  
  err =rep(0,K)
  
  err_true =rep(0,K)
  
  err[1] =1e10
  
  k = 2
  
  ind = 0
  
  
  trans_sd = 0.5
  
  sigma_X1 = 0.5
  
  while(k<K && ind ==0){
    
    
    
    
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
    
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:n){
          M_i = Mean_X[[t]][i,]
          M_j = Mean_X[[t]][j,]
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
    
    
    mean_beta = rnorm(n=1, mean= mean_beta_new,sd = sigma_beta_new)
    
    ### update of Sigma_X, Mean_X
    inde = 1:n
    
    for (i in inde[1:n]){
      for (t in 1:T){
        for (j in inde[1:n]){
          if (which(j==inde) > which(i==inde)){
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
          } else if (which(j==inde) < which(i==inde))
          {
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
          }
          
          if (j!= i){
            V_X_cumulative[[t]][[i]] = V_X_cumulative[[t]][[i]] + sigma^(-2)*( M_j %*% t(M_j) ) *alpha
            M_X_cumulative[[t]][i,] = M_X_cumulative[[t]][i,] +  sigma^(-2)*(Y[[t]][i,j]-mean_beta) * M_j *alpha
          }
        }
        
        if(t == 1)
        { Sigma_X_new[[t]][[i]] = ginv(sigma_X1^(-2)*diag(d)+rho^2*trans_sd^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (rho*trans_sd ^(-2)*Mean_X[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        
        if(1<t  && t<T)
        { Sigma_X_new[[t]][[i]] = ginv((1+rho^2)*trans_sd ^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (rho*trans_sd ^(-2)*Mean_X[[t-1]][i,]+rho*trans_sd ^(-2)*Mean_X[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        if(t==T)
        { Sigma_X_new[[t]][[i]] = ginv(trans_sd ^(-2)*diag(d)+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (rho*trans_sd ^(-2)*Mean_X[[t-1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        Mean_X[[t]][i,] = rmvnorm(n=1, mean = Mean_X_new[[t]][i,],sigma = Sigma_X_new[[t]][[i]])
        
      }
      
    }
      
      b=0
      b0 = 0
      for (i0 in 1:n){
        for (t0 in 1:(T-1)){
          
          EX = sum(Mean_X[[t0+1]][i0,]^2)+sum(Mean_X[[t0]][i0,]^2)-2*t(Mean_X[[t0]][i0,]) %*% Mean_X[[t0+1]][i0,]
          #  sum(diag(x=Sigma_X_new[[t0]][[i0]]))+  sum(diag(x = Sigma_X_new[[t0+1]][[i0]])) 
          
          b = b+ EX
        }
        b0 = b0+ sum(Mean_X[[1]][i0,]^2)
        #+sum(diag(x=Sigma_X_new[[1]][[i0]]))
      }
      
      ded_free = 1-n*(T-1)*d/2
      
      # if (abs(ded_free)<200){
      #   tau_inv_sq = besselK(x=sqrt(b), nu=ded_free+1)/sqrt(b)/besselK(x=sqrt(b), nu=ded_free)-2*ded_free/b 
      # }  else{
      #   tau_inv_sq =  exp(besselK.nuAsym(x=sqrt(b), nu=-ded_free-1,k.max=0,log=T)-besselK.nuAsym(x=sqrt(b), nu=-ded_free,k.max=5,log=T))/sqrt(b)-2*ded_free/b
      # }
      # 
      

      
      tau_inv_sq =  rgig(n=1,  lambda=ded_free, chi=b, psi=1)
      
      trans_sd = as.numeric(sqrt(tau_inv_sq))
      
      sg_inve_sq = rinvgamma(n=1,shape= (n*d+1)/2, rate=(b0+1)/2)
      
      sigma_X1 = as.numeric(sqrt(sg_inve_sq))
      
      
    pred_mean = rep(T*(n-1)*n,0)
    res = rep(T*(n-1)*n,0)
    res_true  = rep(T*(n-1)*n,0)
    
    r=1
   
 #   cat(trans_sd,'\n')
if(!is.null(X_true))    {
    for (t in 1:T){
      for (i in 1:n){
        for (j in 1:n){
          if (j != i){
            pred_mean[r] = mean_beta_new+t(Mean_X_new[[t]][i,])%*% Mean_X_new[[t]][j,]
            res[r] = Y[[t]][i,j]
            res_true[r] = t(X_true[[t]][i,]) %*% X_true[[t]][j,]
            r=r+1
          }
        }
      }
    }
} else if (is.null(X_true)){
  for (t in 1:T){
    for (i in 1:n){
      for (j in 1:n){
        if (j != i){
          pred_mean[r] = mean_beta_new+t(Mean_X_new[[t]][i,])%*% Mean_X_new[[t]][j,]
          res[r] = Y[[t]][i,j]
          r=r+1
        }
      }
    }
  }
  }
    
    err[k]=sqrt(sum((res-pred_mean)^2/T/n/(n-1)))
    
    err_true[k]=sqrt(sum((res_true-pred_mean)^2/T/n/(n-1)))
    
    #abs(err[k]-err[k-1]) <gap||
    
    
    
    if( k>max_iter){
      ind = 1
      
      cat(err[k],'\n')
      
      cat('iteration:',k-1,'\n')
    }    else{
      
      

      
      if (k > (max_iter/2)){
        for(t in 1:T){
          X_out[[t]] =  X_out[[t]]+Mean_X[[t]]/(max_iter/2)
        }
      

      beta_out = beta_out +mean_beta/(max_iter/2)
      }
      # Mean_X = Mean_X_new
      # 
      # Sigma_X = Sigma_X_new
      # 
      # mean_beta = as.numeric( mean_beta_new)
      # 
      # sigma_beta = sigma_beta_new
      
      
      
      k = k+1
    }
    
  }
  
  return(list(err=err[2:(k-1)],trans_sd=trans_sd, Mean_X= X_out,   mean_beta  = beta_out,  iter = k-1 ))
}
