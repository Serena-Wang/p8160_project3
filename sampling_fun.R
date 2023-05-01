
#' log posterior for sigma squared
#' @param sigma_squared double
#' @param Y full response vector in shape N x 1 where N is the total number of observations in the data set
#' @param D full design matrix for random effects for all hurricanes. Its shape is N x 5.
#' @param B full random effects coefficients matrix in shape 5 x H where H is the total number of hurricanes
#' @param X full design matrix for fixed effects for all hurricanes. Its shape is H x 3. 
#' @param gamma fixed effects coefficients matrix in shape 3 x 1. 
#' 
#' @return log posterior value for sigma squared with the given matrices
sigma_squared_log_posterior <- function(sigma_squared, Y, D, B, X, gamma, H, ith_hurricane_idx){
  
  N = length(Y)
  
  mu_H <- rep(NA, N)
  for (i in 1:H){
    curr_hurricane_idx <- ith_hurricane_idx[[i]]
    mu_i <- D[curr_hurricane_idx,,drop = FALSE] %*% B[,i,drop = FALSE]
    mu_H[curr_hurricane_idx] <- mu_i
  }
  
  mu_H <- as.matrix(mu_H)
  
  Y_minus_mu <- Y-mu_H- X%*%gamma
  
  if (sigma_squared>0){
    sigma <- sqrt(sigma_squared)
    
  } else{
    # avoid -inf :) 
    return(-1e10)
  }
  
  log_likelihood <- (-N)* log(sigma) + (t(Y_minus_mu)%*%(Y_minus_mu)/(-2*sigma_squared))
  log_prior <- -log(sigma*pi*(sigma_squared + 100))
  return (log_likelihood+log_prior)
}


#' calculate next value for sigma squared in metropolis hasting random walk  
#' https://towardsdatascience.com/from-scratch-bayesian-inference-markov-chain-monte-carlo-and-metropolis-hastings-in-python-ef21a29e25a
#' @param curr current sigma squared value
#' @param a window size in transition distribution Unif(curr-a, curr+a)
#' @param Y full response vector in shape N x 1 where N is the total number of observations in the data set
#' @param D full design matrix for random effects for all hurricanes. Its shape is N x 5.
#' @param B full random effects coefficients matrix in shape 5 x H where H is the total number of hurricanes
#' @param X full design matrix for fixed effects for all hurricanes. Its shape is H x 16. 
#' @param gamma fixed effects coefficients matrix in shape 16 x 1. 
#' @param H integer number of hurricanes
#' @param ith_hurricane_idx list of list of integers indexes for observations for the ith hurricane
#' 
#' @return next value for sigma squared in metropolis hasting random walk  
metropolis_hasting_sigma_squared_next_val <- function(curr, a, Y, D, B, X, gamma, H, ith_hurricane_idx){
  
  prop <- curr + (runif(1)-0.5)*2*a
  
  prop_post <- sigma_squared_log_posterior(prop, Y, D, B, X, gamma, H, ith_hurricane_idx)
  
  curr_post <- sigma_squared_log_posterior(curr, Y, D, B, X, gamma, H, ith_hurricane_idx)
  
  posterior_ratio <- exp(prop_post - curr_post)
  
  # acceptance prob
  alpha = min(posterior_ratio, 1)
  
  # random sample from unif(0,1)
  u <- runif(n = 1)
  
  if (u < alpha) {
    # accept
    return (prop)
  } else {
    return (curr)
  }
}


metropolis_hasting_sigma_squared <- function(start_sigma_squared, 
                                             numtocompute, 
                                             a = 2,
                                             Y, D, B, X, gamma, H, ith_hurricane_idx){
  chain <- rep(NA, numtocompute)
  chain[1] <- start_sigma_squared
  for (i in 2:numtocompute) {
    #print(i)
    chain[i] <- metropolis_hasting_sigma_squared_next_val(chain[i-1],a,Y, D, B, X, gamma, H, ith_hurricane_idx)
  }
  return (chain)
}


# to get B transpose --> need to transpose B
# B is 5 x H

#' calculate conditional posterior for B
#' 
#' @param D full design matrix for random effects for all hurricanes in shape N x 5.
#' @param sigma_squared a double for sigma squared value
#' @param A inverse of variance of beta's in shape 5 x 5.
#' @param X full design matrix for fixed effects for all hurricanes. Its shape is H x 16. 
#' @param gamma fixed effects coefficients matrix in shape 16 x 1. 
#' @param mu  mean of beta's in shape 5 x 1
#' @param H integer number of hurricanes
#' @param ith_hurricane_idx list of list of integers indexes for observations for the ith hurricane
#' 
#' @return posterior of B in shape 5 x H based on the given values
#'  
B_post <- function(D, sigma_squared, A, X, gamma, mu, H, ith_hurricane_idx){
  B <- NULL
  total_obs <- nrow(D)
  for (i in 1:H){
    curr_hurricane_idx <- ith_hurricane_idx[[i]]
    n_i <- length(curr_hurricane_idx)
    X_i <- X[curr_hurricane_idx,,drop = FALSE]
    Y_i <- Y[curr_hurricane_idx,,drop = FALSE]
    # retain matrix form
    D_i <- D[curr_hurricane_idx,,drop = FALSE] 
    
    M_i <- t(D_i)%*% ((1/sigma_squared)*diag(n_i)) %*% D_i + A
    
    N_i<- t(D_i)%*% ((1/sigma_squared)*diag(n_i)) %*% Y_i - t(D_i)%*% (1/sigma_squared*diag(n_i)) %*% X_i%*%gamma + A%*%mu
    
    B_i <- mvrnorm(1, mu = solve(M_i) %*% N_i, Sigma = solve(M_i))
    B <- cbind(B, B_i)
  }
  
  return (B)
}


#' calculate conditional posterior for mu
#'
#' @param A inverse of variance of beta's in shape 5 x 5.
#' @param V variance-covariance matrix for mu in shape 5 x 5
#' @param B full random effects coefficients matrix in shape 5 x H where H is the total number of hurricanes
#' 
#' @return posterior of mu in shape 5 x 1 based on the given values
mu_post <- function(A, V, B){
  M <- H* A - solve(V)
  N <- 0
  for (i in 1:H){
    N <- N + A %*% B[,i,drop = FALSE]
  }
  mu <-  mvrnorm(1, mu = solve(M) %*% N, Sigma = solve(M))
  
  return(mu)
}


#' calculate conditional posterior for A (inverse of Sigma)
#'
#' @param nu integer for degree of freedom
#' @param B full random effects coefficients matrix in shape 5 x H where H is the to
#' @param mu  mean of beta's in shape 5 x 1
#' @param S scale matrix for A in shape 5 x 5.
#' @param H integer number of hurricanes
#' 
#' @return posterior of A in shape 5 x 5 based on the given values
A_post <- function(nu, B, mu, S, H){
  df <- H + nu
  
  sum = 0
  for (i in 1:H){
    sum = sum + (B[,i,drop = FALSE] - mu) %*% t((B[,i,drop = FALSE] - mu))
  }
  
  scale <- solve(S+ sum)
  
  A <- matrixsampling::rwishart(n = 1, nu = df, Sigma = scale, checkSymmetry = F)
  
  A <- array(A[, , 1], dim=c(5, 5)) 
  return (A)
}


#' calculate conditional posterior for gamma
#' 
#' @param X full design matrix for fixed effects for all hurricanes. Its shape is H x 16. 
#' @param Y full response vector in shape N x 1 where N is the total number of observations in the data set
#' @param D full design matrix for random effects for all hurricanes in shape N x 5.
#' @param B full random effects coefficients matrix in shape 5 x H where H is the total number of hurricanes
#' @param sigma_squared  a double for sigma squared value
#' @param H integer number of hurricanes
#' @param ith_hurricane_idx list of list of integers indexes for observations for the ith hurricane
#' 
#' @return  posterior of gamma in shape 16 x 1 based on the given values
gamma_post <- function(X,Y, D,B, sigma_squared, H, ith_hurricane_idx){
  M <- 0
  N <- 0
  for (i in 1:H){
    curr_hurricane_idx <- ith_hurricane_idx[[i]]
    n_i <- length(curr_hurricane_idx)
    X_i <- X[curr_hurricane_idx,,drop = FALSE]
    Y_i <- Y[curr_hurricane_idx,,drop = FALSE]
    D_i <- D[curr_hurricane_idx,,drop = FALSE]
    M = M + t(X_i) %*% ((1/sigma_squared)*diag(n_i)) %*% X_i
    
    N = N + t(X_i) %*% ((1/sigma_squared)*diag(n_i))  %*% Y_i - t(X_i) %*% (1/sigma_squared*diag(n_i)) %*% D_i %*% B[,i,drop = FALSE]
    
  }
  M = M+ 400*diag((dim(X)[2]))
  
  gamma <- mvrnorm(1, mu = solve(M) %*% N, Sigma = solve(M)) %>% as.matrix()
  
  return(gamma)
}

test_sigma_squared_post <- function(Y, D, B , X, gamma, H, ith_hurricane_idx){
  
  shape <- (nrow(D)-1)/2
  
  mu_H <- rep(NA, N)
  for (i in 1:H){
    curr_hurricane_idx <- ith_hurricane_idx[[i]]
    mu_i <- D[curr_hurricane_idx,,drop = FALSE] %*% B[,i,drop = FALSE]
    mu_H[curr_hurricane_idx] <- mu_i
  }
  
  mu_H <- as.matrix(mu_H)
  
  Y_minus_mu <- Y-mu_H- X%*%gamma
  
  rate <-t(Y_minus_mu)%*%(Y_minus_mu)/2
  
  
  sigma_squared <- rinvgamma(1, shape = shape, scale = rate)
  
  return(sigma_squared)
}



hurricane_gibbs_sampling <- function(start_sigma_squared, 
                                     start_A,
                                     start_gamma, 
                                     start_mu,
                                     V,
                                     S,
                                     X,
                                     D,
                                     Y,
                                     nitr = 100,
                                     nu = 5, 
                                     H, 
                                     ith_hurricane_idx){
  B_post_results <- list()
  mu_post_reulsts <- list()
  A_post_results <- list()
  gamma_post_results <- list()
  sigma_squared_post_results <- list()
  chains<-list()
  
  # B first
  B <- B_post(D, start_sigma_squared, start_A, X, start_gamma, start_mu, H, ith_hurricane_idx)
  mu <- mu_post(A = start_A, V= V, B = B)
  A <- A_post(nu = nu, B = B, mu = mu, S = S, H = H)
  gamma <- gamma_post(X = X, 
                      Y = Y, 
                      D = D,
                      B = B,
                      sigma_squared = start_sigma_squared,
                      H = H, 
                      ith_hurricane_idx = ith_hurricane_idx)
  
  sigma_squared_chain <- metropolis_hasting_sigma_squared(
    start_sigma_squared = start_sigma_squared, 
    numtocompute = 500, 
    a = 0.8, 
    Y = Y, 
    D = D, 
    B = B, 
    X = X, 
    gamma = gamma,
    H = H, 
    ith_hurricane_idx = ith_hurricane_idx)
  
  sigma_squared <- mean(sigma_squared_chain[401:500])
  
  B_post_results[[1]] <- B
  mu_post_reulsts[[1]] <- mu
  A_post_results[[1]] <- A
  gamma_post_results[[1]] <-gamma
  sigma_squared_post_results[[1]] <-sigma_squared
  chains[[1]] <- sigma_squared_chain
  
  for (i in 2:nitr){
    
    print(i)
    
    B <- B_post(D = D, 
                sigma_squared = sigma_squared_post_results[[i-1]], 
                A = A_post_results[[i-1]], 
                X = X, 
                gamma = gamma_post_results[[i-1]], 
                mu = mu_post_reulsts[[i-1]],
                H = H, 
                ith_hurricane_idx = ith_hurricane_idx)
    
    mu <- mu_post(A = A_post_results[[i-1]], 
                  V= V, 
                  B = B)
    A <- A_post(nu = nu, 
                B = B, 
                mu = mu, 
                S = S,
                H = H)
    
    gamma <- gamma_post(X = X,
                        Y = Y, 
                        D = D,
                        B = B,
                        sigma_squared = sigma_squared_post_results[[i-1]],
                        H = H, 
                        ith_hurricane_idx = ith_hurricane_idx)
    
    # sigma_squared on its own random walk
    sigma_squared_chain <- metropolis_hasting_sigma_squared(
      start_sigma_squared = sigma_squared_post_results[[i-1]], 
      numtocompute = 500, 
      a = 0.8, 
      Y = Y, 
      D = D, 
      B = B, 
      X = X, 
      gamma = gamma,
      H = H, 
      ith_hurricane_idx = ith_hurricane_idx)
    
    sigma_squared <- mean(sigma_squared_chain[401:500])
    
    # collect results
    B_post_results[[i]] <- B
    mu_post_reulsts[[i]] <- mu
    A_post_results[[i]] <- A
    gamma_post_results[[i]] <-gamma
    sigma_squared_post_results[[i]] <-sigma_squared
    chains[[i]] <- sigma_squared_chain
  }
  
  return(list(
    B = B_post_results,
    mu = mu_post_reulsts,
    A = A_post_results,
    gamma = gamma_post_results,
    sigma_squared = sigma_squared_post_results,
    chains = chains
  ))
  
}

