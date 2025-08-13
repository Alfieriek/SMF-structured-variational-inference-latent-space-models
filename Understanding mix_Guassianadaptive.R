source('Understanding simu Guassian Simple MEAN FIELD CASE.R')
# From simu_Gaussian_simple
Y # List of Socio Matrix
True_X <- X # List of True Latent Variables
True_beta <- 0 # True value for parameter beta

rm(list = setdiff(ls(), c("Y","True_X","True_beta"))) # Removing variables not needed


# parameters of beta prior
mean_beta_prior = 0    #prior mean of beta, 
sigma_beta_prior = 10   #prior sd of beta, 


# sigma: sd for the likelihood
sigma = 0.1 #sd for likelihood of y





#----------------------------------------------------------------------------
# Gaussian network via SMF variational inference
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





# mix_Gaussian_adaptive = function(Y ,mean_beta_prior, d=2, sigma_beta_prior,rho=1, sigma, X_true=NULL, gap = 1e-3, max_iter=10, alpha=0.95,global_prior='Gamma', IG=1, IG_2 =1/2){
#   
#   
#   
Y 
mean_beta_prior
d=2
sigma_beta_prior
rho=1
sigma
X_true=NULL
gap = 1e-3
max_iter=10
alpha=0.95
global_prior='Gamma'
IG=1
IG_2 =1/2

  
  
  
#----------  
# Extracting the Number of timesteps "T" and dimension of socio matrices "n"
#----------
  
  T = length(Y)
  
  n = length(Y[[1]][1,])
  
  
  #----------  
  # Initializing Variational Parameters
  #----------  
  
  #### target parameters
  
  mean_beta = 0;    #mean of beta, 
  
  sigma_beta = sqrt(2);  #sd of beta, 
  
  Mean_X = vector("list", T); #Mean of X, Mean_X[[t]][i,] 
  
  
  Sigma_X = replicate(n=T, expr=list()); # covariance matrix of X, Sigma_x[[t]][[i]] is the covariance matrix Sigma_{it}
  
  
  for (t in 1:T){
    Mean_X[[t]] = matrix(rep(0,n*d),nrow = n)
    Sigma_X[[t]] = vector("list", n)
    for (i in 1:n){
      Mean_X[[t]][i,] =  rmvnorm(n=1,mean=rep(0,d), sigma=0.1*diag(d)) # randomly initialization for mu_{it}
      Sigma_X[[t]][[i]] =  diag(d)          # initialization for Sigma_{it}
    }
  }
  
  
  #--------------------------------Algorithm---------------------
  K= 500 # What is this? Number of iterations? # "max_iter" is not used!!!. 
  
  err =rep(0,K) # Error vector to compare with "gap"
  
  err_true =rep(0,K) # True Error vector ?
  
  
  err[1] =1e10 # Replacing First Error value to start loop
  
  k = 2 # Second iterate as initialized value is k = 1
  
  ind = 0 # index/indices? 
  
  #----------  
  # Updated Variational Parameters variables
  #----------  
  
  mean_beta_new = 0;    
  
  sigma_beta_new = 0  
  
  Mean_X_new = vector("list", T) # List
  
  Sigma_X_new = replicate(n=T, expr=list()) # List of Lists
  
  Cross_cov = replicate(n=T, expr=list()) # List of Lists
  
  for (t in 1:T){
    Mean_X_new[[t]] = matrix(rep(0,n*d),nrow = n)
    Sigma_X_new[[t]] = vector("list", n)
    Cross_cov[[t]] = vector("list", n)
    for (i in 1:n){
      Sigma_X_new[[t]][[i]] =  0*diag(d)
      Cross_cov[[t]][[i]] = 0*diag(2*d)
    }
  }
  
  #----------  
  # 
  #----------  
  
  trans_sd = 0.5 # Prior for transitional state
  
  sigma_X1 = 0.5 # Prior of initial state
  
  #-- Redundant
  Mean_X_new = vector("list", T) # List of List
  
  Sigma_X_new = replicate(n=T, expr=list()) # List of Lists
  #-- Redundant
  
  Forward.Mean.new = vector("list", T) # List
  
  Backward.Mean.new = vector("list", T) # List
  
  Forward.var.new = replicate(n=T, expr=list()) # List of Lists
  
  Backward.var.new = replicate(n=T, expr=list()) # List of Lists
  
  
  for (t in 1:T){
    Mean_X_new[[t]] = matrix(rep(0,n*d),nrow = n)
    Forward.Mean.new[[t]] = matrix(rep(0,n*d),nrow = n)
    Backward.Mean.new[[t]] = matrix(rep(0,n*d),nrow = n)
    Sigma_X_new[[t]] = vector("list", n)
    Forward.var.new[[t]] = vector("list", n)
    Backward.var.new[[t]] = vector("list", n)
    for (i in 1:n){
      Sigma_X_new[[t]][[i]] =  0*diag(d)       
      Forward.var.new[[t]][[i]] =  0*diag(d)
      Backward.var.new[[t]][[i]] =  0*diag(d)
    }
  }
  
  
  while(k<K && ind ==0){
    V_beta_cumulative = 0 # Not USED!!!!
    
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
    
    # Calculate auxiliary values
    
    #--- Page 49.  This is the double sum part without sigma^{-2} 
    for (t in 1:T){
      for( i in 1:n){
        for (j in 1:n){
          M_i = Mean_X[[t]][i,]
          M_j = Mean_X[[t]][j,]
          V_i = Sigma_X[[t]][[i]]
          V_j = Sigma_X[[t]][[j]]
          if (j!= i){ 
            M_beta_cumulative = M_beta_cumulative + (Y[[t]][i,j]- t(M_i) %*% M_j)*alpha 
          }
        }
      }
    }
    
    
    
    ### update of sigma_beta
    
    sigma_beta_new = (sigma_beta_prior^(-2)+ T*n*(n-1)*sigma^(-2) )^(-1/2) # Page 49, Where is alpha??
    
    
    ### update of mean_beta
    
    mean_beta_new = as.numeric(sigma_beta_new^2* (sigma_beta_prior^(-2)*mean_beta_prior + sigma^(-2)*M_beta_cumulative) )
    
    
    
    ### update of Sigma_X, Mean_X
    
    
    
    for (i in 1:n){
      
      
      #----------------------------------------------------------------------- Computing Cumulative Potentials Mean and Var.
      
      for (t in 1:T){
        
        for (j in 1:n){
          
          
          
          if (j>i){
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X[[t]][j,]
            V_i = Sigma_X[[t]][[i]]
            V_j = Sigma_X[[t]][[j]]
          } else if (j<i)
          {
            M_i = Mean_X[[t]][i,]
            M_j = Mean_X_new[[t]][j,]
            V_i = Sigma_X[[t]][[i]]
            V_j = Sigma_X_new[[t]][[j]] 
          }
          
          if (j!= i){
            V_X_cumulative[[t]][[i]] = V_X_cumulative[[t]][[i]] +  sigma^(-2)*(M_j %*% t(M_j) + V_j) *alpha
            M_X_cumulative[[t]][i,] = M_X_cumulative[[t]][i,] +  sigma^(-2)*(Y[[t]][i,j]-mean_beta_new) * M_j *alpha
          }
        }
      }
      
      #-----------------------------------------------------------------------
      # Computing messages Forward.mean.new = \mu_{i, t \to t+1 }  and forward.var.new = \Sigma_{i, t \to t+1 } 
      Forward.var.new[[1]][[i]] = -rho^(-2)*trans_sd^(4)*(sigma_X1^(-2)*diag(d)+rho^2*trans_sd^(-2)*diag(d)+V_X_cumulative[[1]][[i]])
      
      Forward.Mean.new[[1]][i,] = -rho^(-1)*trans_sd^2*
        (M_X_cumulative[[1]][i,])
      

      Backward.var.new[[T]][[i]] = -rho^(-2)*trans_sd^(4)*(trans_sd^(-2)*diag(d)+V_X_cumulative[[T]][[i]])
      
      Backward.Mean.new[[T]][i,] = -rho^(-1)*trans_sd^2*
        (M_X_cumulative[[T]][i,])
      
      
      for (t2 in 2:(T-1)){

        Forward.var.new[[t2]][[i]] = -rho^(-2)*trans_sd^(4)*(ginv(Forward.var.new[[t2-1]][[i]])+(1+rho^2)*trans_sd^(-2)*diag(d)+V_X_cumulative[[t2]][[i]])
        
        Forward.Mean.new[[t2]][i,] = -rho^(-1)*trans_sd^2*(ginv(Forward.var.new[[t2-1]][[i]])%*%  Forward.Mean.new[[t2-1]][i,]  +M_X_cumulative[[t2]][i,])
        
        
        Backward.var.new[[T+1-t2]][[i]] = -rho^(-2)*trans_sd^(4)*(ginv(Backward.var.new[[T+2-t2]][[i]])+(1+rho^2)*trans_sd^(-2)*diag(d)+V_X_cumulative[[T+1-t2]][[i]])
        
        Backward.Mean.new[[T+1-t2]][i,] = -rho^(-1)*trans_sd^2
        (ginv(Backward.var.new[[T+2-t2]][[i]])%*% Backward.Mean.new[[T+2-t2]][i,] +M_X_cumulative[[T+1-t2]][i,])
      }
      
      #-----------------------------------------------------------------------
      
      # MYBW UPDATES PG 14
      #----------------------------------------------------------------------- This is the first equation from eq(20)
      for(t in 1:T){
        if(t == 1)
        { Sigma_X_new[[t]][[i]] = ginv(sigma_X1^(-2)*diag(d)+rho^2*trans_sd^(-2)*diag(d)+ginv(Backward.var.new[[t+1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Backward.var.new[[t+1]][[i]])%*% Backward.Mean.new[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        if(t==T)
        { Sigma_X_new[[t]][[i]] = ginv(trans_sd^(-2)*diag(d)+ginv(Forward.var.new[[t-1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Forward.var.new[[t-1]][[i]])%*%  Forward.Mean.new[[t-1]][i,]+M_X_cumulative[[t]][i,])
        }
        
        if(1<t  && t<T)
        { Sigma_X_new[[t]][[i]] = ginv((1+rho^2)*trans_sd^(-2)*diag(d)+ginv(Forward.var.new[[t-1]][[i]])+ginv(Backward.var.new[[t+1]][[i]])+V_X_cumulative[[t]][[i]])
        
        Sigma_X_new[[t]][[i]] = (Sigma_X_new[[t]][[i]] +t(Sigma_X_new[[t]][[i]]))/2
        
        Mean_X_new[[t]][i,] = Sigma_X_new[[t]][[i]]%*%
          (ginv(Forward.var.new[[t-1]][[i]])%*%  Forward.Mean.new[[t-1]][i,]+
             ginv(Backward.var.new[[t+1]][[i]])%*% Backward.Mean.new[[t+1]][i,]+M_X_cumulative[[t]][i,])
        }
        
      }
      #----------------------------------------------------------------------- 
      
      #----------------------------------------------------------------------- This is the second equation from eq(22)
      Cross_cov[[1]][[i]][1:d,1:d] = sigma_X1^(-2)*diag(d)+trans_sd^(-2)*diag(d)+V_X_cumulative[[1]][[i]]
      Cross_cov[[1]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[1]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[1]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+ginv(Backward.var.new[[3]][[i]])+(trans_sd^(-2))*diag(d)+V_X_cumulative[[2]][[i]]
      Cross_cov[[1]][[i]] = ginv(Cross_cov[[1]][[i]])
      
      for ( t in 2:(T-2)){
        Cross_cov[[t]][[i]][1:d,1:d] = (trans_sd^(-2))*diag(d)+ ginv(Forward.var.new[[t-1]][[i]]) +(trans_sd^(-2))*diag(d)+V_X_cumulative[[t]][[i]]
        Cross_cov[[t]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
        Cross_cov[[t]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
        Cross_cov[[t]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+ginv(Backward.var.new[[t+2]][[i]])+(trans_sd^(-2))*diag(d)+V_X_cumulative[[t+1]][[i]]
        Cross_cov[[t]][[i]] = ginv(Cross_cov[[t]][[i]])
      }
      
      Cross_cov[[T-1]][[i]][1:d,1:d] = (trans_sd^(-2))*diag(d)+ ginv(Forward.var.new[[T-2]][[i]]) +(trans_sd^(-2))*diag(d)+V_X_cumulative[[T-1]][[i]]
      Cross_cov[[T-1]][[i]][1:d,(d+1):(2*d)] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[T-1]][[i]][(d+1):(2*d),1:d] = -(trans_sd^(-2))*diag(d)
      Cross_cov[[T-1]][[i]][(d+1):(2*d),(d+1):(2*d)] = (trans_sd^(-2))*diag(d)+V_X_cumulative[[T]][[i]]
      Cross_cov[[T-1]][[i]] = ginv(Cross_cov[[T-1]][[i]])
      #----------------------------------------------------------------------- 
      
      
      
      
      
      
      b=0
      b0 = 0
      for (i0 in 1:n){
        for (t in 1:(T-1)){
          
          EX = sum(Mean_X_new[[t+1]][i0,]^2)+sum(Mean_X_new[[t]][i0,]^2)-2*t(Mean_X_new[[t]][i0,]) %*% Mean_X_new[[t+1]][i0,]+
            sum(diag(x= as.matrix(Cross_cov[[t]][[i0]][1:d,1:d]+ Cross_cov[[t]][[i0]][(d+1):(2*d),(d+1):(2*d)] -
                                    2*Cross_cov[[t]][[i0]][1:d,(d+1):(2*d)])))
          
          # EX = sum(Mean_X_new[[t+1]][i,]^2)+sum(Mean_X_new[[t]][i,]^2)-2*t(Mean_X_new[[t]][i,]) %*% Mean_X_new[[t+1]][i,]+
          #   sum(diag(x=Sigma_X_new[[t]][[i]]))+  sum(diag(x = Sigma_X_new[[t+1]][[i]])) 
          
          b = b+ EX
        }
        b0 = b0 + sum(Mean_X_new[[1]][i0,]^2)+sum(diag(x=Sigma_X_new[[1]][[i0]]))
      }
      
      # ded_free = 1-n*(T-1)*d/2
      # 
      # if (abs(ded_free)<200){
      #   tau_inv_sq = besselK(x=sqrt(b), nu=ded_free+1)/sqrt(b)/besselK(x=sqrt(b), nu=ded_free)-2*ded_free/b 
      # }  else{
      #   tau_inv_sq =  exp(besselK.nuAsym(x=sqrt(b), nu=-ded_free-1,k.max=5,log=T)-besselK.nuAsym(x=sqrt(b), nu=-ded_free,k.max=5,log=T))/sqrt(b)-2*ded_free/b
      # }
      
      if (global_prior =='Gamma'){
        
        gig_deg = -(T-1)*d*n/2+1/2-IG/2
        
        tau_inv_sq = sqrt(2*IG_2)*exp(besselK.nuAsym(x=sqrt(2*IG_2*b), nu=-gig_deg-1,k.max=5,log=T)-besselK.nuAsym(x=sqrt(2*IG_2*b), nu=-gig_deg,k.max=5,log=T))/sqrt(b)-2*gig_deg/b
        
      } else if (global_prior =='Cauthy'){
        v0 = 1/(1+1/trans_sd^2)
        
        tau_inv_sq = (d*(T-1)*n/2+1/2)/(v0+b/2)
        
      }
      
      trans_sd = as.numeric(1/sqrt(tau_inv_sq))
      
      #  cat(trans_sd,'\n')
      
      sg_inve_sq = n*d+1/(b0+1)
      
      sigma_X1 = as.numeric(1/sqrt(sg_inve_sq))
    }
    
    
    
    
    
    pred_mean = rep(T*(n-1)*n,0)
    res = rep(T*(n-1)*n,0)
    res_true  = rep(T*(n-1)*n,0)
    res_Y = rep(T*(n-1)*n,0)
    
    r=1
    
    if(!is.null(X_true))    {
      for (t in 1:T){
        for (i in 1:n){
          for (j in 1:n){
            if (j != i){
              pred_mean[r] = mean_beta_new+t(Mean_X_new[[t]][i,])%*% Mean_X_new[[t]][j,]
              res[r] = mean_beta+t(Mean_X[[t]][i,])%*% Mean_X[[t]][j,]
              res_Y[r] = Y[[t]][i,j]
              res_true[r] = t(X_true[[t]][i,]) %*% X_true[[t]][j,]+beta
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
              res[r] = mean_beta+t(Mean_X[[t]][i,])%*% Mean_X[[t]][j,]
              res_Y[r] = Y[[t]][i,j]
              r=r+1
            }
          }
        }
      }
    }
    
    err[k]=sqrt(sum((res-pred_mean)^2/T/n/(n-1)))
    
    err_true[k]=sqrt(sum((res_true-pred_mean)^2/T/n/(n-1)))
    
    if((abs(err[k]) <gap ) || k>50){
      ind = 1
      cat(err[k],'\n')
      
      
      cat('iteration:',k-1,'\n')
    }    else{
      #err[k]-err[k-1] >-0.001
      
      #err[k]-err[k-1]
      
      Mean_X = Mean_X_new
      
      Sigma_X = Sigma_X_new
      
      mean_beta = as.numeric( mean_beta_new)
      
      sigma_beta = sigma_beta_new
      
      
      
      
      
      k = k+1
    }
    
  }
  
#   return(list(err=err[2:(k-1)],trans_sd=trans_sd, Mean_X= Mean_X, Sigma_X= Sigma_X,  mean_beta  = mean_beta, sigma_beta =sigma_beta , iter = k-1 ))
# }

