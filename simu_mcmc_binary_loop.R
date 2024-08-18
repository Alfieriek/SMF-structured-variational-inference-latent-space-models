rm(list = ls())
set.seed(2023)
library(MASS)
library(mvtnorm)
library(igraph)
library(Bessel)

file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(file_location)
source('helper.R')
source('mean_field_DN_adaptive.R')
source('mix_DN_adaptive.R')
source('mcmc_DN_adaptive.R')

# Simulation Function
run_simulation <- function(simulation_count = 25, n_values = c(10, 20, 50), tau_values = c(0.01,0.05,0.1)) {
  # Dataframe to hold results

  
  # Loop over different n and tau values
  for(n in n_values) {
    for(tau in tau_values) {
      # Simulation loop
      for(sim_id in 1:simulation_count) {
        # ...
        results <- data.frame(
          sim_id = integer(),
          method = character(),
          iterations = integer(),
          n = integer(),
          tau = numeric(),
          time = numeric(),
          error = numeric()
        )

        T = 100   # T = 100;
        
        rho = 0.5
        
        beta = -1 # 
        
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
        
        sigma_beta_prior = sqrt(10)   #prior sd of beta_0, 
        
        
        gap = 1e-4
        
        min_iter = 1
        # Loop over different iteration counts
        for(max_iter in c(50, 100, 200)) {
          # Mean Field Model Fitting
          start_time <- Sys.time()
          MF_list_adaptive =mcmc_DN_adaptive (Y = Y,d=d,mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                                     gap = gap,max_iter=max_iter)
          end_time <- Sys.time()
          
          # Compute prediction error
          pred_mean_MF_adaptive = matrix(rep(0,T*n*(n-1)),nrow = T)
          MF_adaptive_sm = rep(0,T)
          res =  matrix(rep(0,T*n*(n-1)),nrow = T)
          
          for (t in 1:T){
            r=1
            for (i in 1:n){
              for (j in 1:n){
                if (j!=i){
                 # Compute squared errors for new MF cases
                  pred_mean_MF_adaptive[t,r] = t(MF_list_adaptive$Mean_X[[t]][i,])%*% MF_list_adaptive$Mean_X[[t]][j,]+MF_list_adaptive$mean_beta
                  res[t,r] = t(X[[t]][i,]) %*% X[[t]][j,]+beta
                  r=r+1
                }
              }
            }
            MF_adaptive_sm[t] = cor(res[t,], pred_mean_MF_adaptive[t,],method = "pearson")
          }
          
          error <- mean(MF_adaptive_sm)
          
          # Append results
          results <- rbind(results, data.frame(
            sim_id = sim_id,
            method = "MF",
            iterations = max_iter,
            n = n,
            tau = tau,
            time = as.numeric(end_time - start_time, units = "secs"),
            error = error
          ))
          
        }
          
          # Structured Mean Field Model Fitting
          start_time <- Sys.time()
          Mix_list_adaptive = mix_DN_adaptive(Y = Y, d=d, mean_beta_prior = mean_beta_prior, sigma_beta_prior = sigma_beta_prior,
                                              gap = gap)
          end_time <- Sys.time()
          
          # Compute prediction error
          pred_mean_Mix_adaptive = matrix(rep(0,T*n*(n-1)),nrow = T)
          Mix_adaptive_sm = rep(0,T)
          res =  matrix(rep(0,T*n*(n-1)),nrow = T)
          
          for (t in 1:T){
            r=1
            for (i in 1:n){
              for (j in 1:n){
                if (j!=i){
                  # Compute squared errors for new MF cases
                  pred_mean_Mix_adaptive[t,r] = t(Mix_list_adaptive$Mean_X[[t]][i,])%*% Mix_list_adaptive$Mean_X[[t]][j,]+Mix_list_adaptive$mean_beta
                  res[t,r] = t(X[[t]][i,]) %*% X[[t]][j,]+beta
                  r=r+1
                }
              }
            }
            Mix_adaptive_sm[t] = cor(res[t,],pred_mean_Mix_adaptive[t,],method = "pearson")
          }
          
          error <- mean(Mix_adaptive_sm)
          
          # Append results
          results <- rbind(results, data.frame(
            sim_id = sim_id,
            method = "SMF",
            iterations = max_iter,
            n = n,
            tau = tau,
            time = as.numeric(end_time - start_time, units = "secs"),
            error = error
          )
          )
          write.table(results, file = "simulation_results_DN.csv", append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
        
      }
    }
  }
}


# Initialize CSV file with column names
write.table(data.frame(
  sim_id = integer(),
  method = character(),
  iterations = integer(),
  n = integer(),
  tau = numeric(),
  time = numeric(),
  error = numeric()
), file = "simulation_results_DN.csv", append = FALSE, sep = ",", col.names = TRUE, row.names = FALSE)

# Run simulation
run_simulation(simulation_count = 25)


