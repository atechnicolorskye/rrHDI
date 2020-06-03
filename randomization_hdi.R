library(glmnet)
library(hdi)
library(mvtnorm)
library(RPtests)

Rcpp::sourceCpp("cRAMSES_V.cpp")
source("RAMSES_V.R") # Randomization Inferenc


get_perm <- function(p){
  ind <- sample(p)
  while(any(ind == 1:p)){
    ind <- sample(p)
  }
  return(diag(length(ind))[ind,])
}


rr_dantzig = function(y, X, n_g, M, a_ind, a_val, n_ind, n_val, g_design, scale_res) {
  # M = sparse M^T obtained via Dantzig selector
  # Test H0_j: beta_j = val_j individually
  
  n = nrow(X)
  p = ncol(X)
  stopifnot(length(y)==nrow(X))
  
  test_ind = c(a_ind, n_ind)
  ac = length(a_ind)
  test_val = c(a_val, n_val)
  na = length(n_val)
  
  t = ac + na
  Y = matrix(y, nrow=n, ncol=t, byrow=F) # n x p
  # Construct binary indicator for H0_j
  A = matrix(0, nrow=t, ncol=p)          # s x p
  for(i in 1:t) { A[i, test_ind[i]] = 1 }
  
  # Run Sqrt LASSO + correction
  sqrt_lasso = RPtests::sqrt_lasso(X, c(y), output_all = TRUE)
  eps = (y - X %*% sqrt_lasso$beta)
  beta_dlasso = sqrt_lasso$beta + 1 / n * M %*% t(X) %*% eps
  if (scale_res == TRUE){
    # Inflates epsilons according to eq. 24 of Zhang and Guang
    eps = eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }
  
  An = (1/sqrt(n)) * A %*% M %*% t(X)   # (s + # inactive) x n.
  Tobs = t(sqrt(n) * (A %*% beta_dlasso - test_val))
  
  Tvals = matrix(0, nrow=0, ncol=t)     # R x (s + # inactive)
  
  R = 1
  while(R < n_g){
    if(g_design=="perm"){
      G = get_perm(n)
    } else if(g_design=="sign"){
      G = diag(sample(c(rep(1, n/2), rep(-1, n/2))))
    }
    if(max(abs(Matrix(G - diag(rep(1,n)))))>0){
      Tvals = rbind(Tvals, t(An %*% G %*% eps))
      R = R +1
    }
  }
  
  get_q = sapply(1:t, function(j){
    quantile(Tvals[,j], c(.025, .975)) /sqrt(n)
  })
  
  ci_a = cbind(beta_dlasso[a_ind]-get_q[2,1:ac], beta_dlasso[a_ind]-get_q[1,1:ac])
  ci_n = cbind(beta_dlasso[n_ind]-get_q[2,(ac+1):t], beta_dlasso[n_ind]-get_q[1,(ac+1):t])
  
  returnList <- list("ci_a" = ci_a,"ci_n" = ci_n)
  return(returnList)
}


