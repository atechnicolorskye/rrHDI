## Reproduce ridge projection results for HDI (Bernoulli 2013)
rm(list=ls())
library("doParallel")
library("fastclime")
library("glmnet")
library("mvtnorm")

library("HDCI") # Buhlmann
source("lasso_inference.r") # Montanari
source("randomization_hdi.R") # Us

registerDoParallel(2)

set.seed(2020)
## DGP
#' (n, p) = (30, 60)   Rand > BLPR, DeLASSO
#'        = (150, 20) or (150, 200) and s =large then BLPR breaks.

# samples
n <- 50 
# dim
p <- 200 # BLPR breaks at 20.
# sparsity
s <- 5  # Set >= 50 and BLPR breaks.

rho = 0.5 
# Sig = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
# S  = |i-j|^r + A*diag(p)
Sig = 4 * diag(p) + rho ^ abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
# Sig = toeplitz(.8^(0:(p-1)))
nsim = 3
# b = 25
Cover = matrix(0,  nrow=0, ncol=3)
Len = matrix(0, nrow=0, ncol=3)

colnames(Cover) = c("BLPR", "DeLASSO", "Randomiz.")
colnames(Len) = c("BLPR", "DeLASSO", "Randomiz.")

# intlen = c()

# for (i in 1:nsim) {
  i = 1
  print(sprintf("--===-  %d/%d iter ==--===",i, nsim))
  beta = rep(0, p)
  beta[1:s] = runif(s, 1/3, 1)
  # shuffle.
  beta = sample(beta)
  # focus on 15 largest
  ind = rev(order(beta))[1:15] # focus on those tests.
  beta_0 = beta[ind]
  
  x = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sig)  # normal covariates.
  ypred = x %*% beta
  
  # sigma = sigma * exp(x %*% beta) / (1 + exp(x %*% beta))
  # sigma = sigma * abs(rowMeans(x^2))
  # y <- x %*% beta + rnorm(n, 0, sigma**2)
  Id = 1:n
  sigma = sqrt(abs(ypred)) * exp(0.05*abs(ypred))
  
  err = sigma * rt(n, df=5)
  plot(ypred, ypred, col="red", type="l")
  points(ypred, ypred + err, col="blue", pch=20)
  y = ypred + err
  
  
  ## Randomization
  print("> Residual Randomization")
  
  # DantzigM
  S <- (t(x) %*% x) / n
  lambda_min <- min(.1, sqrt(log(p) / n))
  n_solve <- 50
  dantzig_M <- fastclime(S, lambda.min=lambda_min, nlambda=n_solve)
  # M = dantzig_M$sigmahat
  M = dantzig_M$icovlist[[length(dantzig_M$icovlist)]]
  S_test = ind
  lam0_vals = beta_0
  
  n = nrow(x)
  p = ncol(x)
  stopifnot(length(y)==nrow(x))
  
  s = length(S_test)
  Y = matrix(y, nrow=n, ncol=s, byrow=F) # make this nxp matrix
  a = matrix(0, nrow=s, ncol=p)  # s x p   binary indicator of  H0_j
  for(i in 1:s) { a[i, S_test[i]] = 1 }
  
  t0 = proc.time()[3]
  # test H0_j: beta_j = v0, for all j in S_test.
  # 1. Create restricted Lasso and debiased Lasso solutions.
  
  beta_lasso = RPtests::sqrt_lasso(x, c(y))
  Eps = (y - x %*% beta_lasso)
  beta_dlasso = beta_lasso + 1 / n * M %*% t(x) %*% Eps
  
  # 2. restricted errors.
  An = (1/sqrt(n)) * a %*% M %*% t(x)   # s x n.
  Tobs = t(sqrt(n) * (a %*% beta_dlasso - lam0_vals))
  
  Tvals = matrix(0, nrow=0, ncol=s)   # R x s
  R = 1
  while(R < 1000){
    G = diag(sample(c(-1, 1), size=n, replace=T))
    if(max(abs(Matrix(G - diag(rep(1,n)))))>0){
      Tvals = rbind(Tvals, t(An %*% G %*% Eps))
      R = R +1
    }
  }
  
  # 3. PVals.
  Pvals = sapply(1:s, function(j) {
    out_pval(list(tobs=Tobs[j], tvals=Tvals[,j]), ret_pval = T)
  })
  
  # 4. CI
  get_q = sapply(1:s, function(j){
    quantile(Tvals[,j], c(.025, .975)) /sqrt(n)
  })
  
  ci = cbind(beta_dlasso[S_test] -get_q[2,], beta_dlasso[S_test]-get_q[1,])

  # out3 = riSingle_hdi_dantzig(y, x, dantzig_M$sigmahat, S_test=ind, lam0_vals=beta_0)
  # cover_ri_M = mean(out3 >=  0.025)
  # len_ri_M = 0 # mean length.
  # 
  # Cover = rbind(Cover, c(cover_blpr, cover_delasso, cover_ri, cover_ri_M))
  # Len = rbind(Len, c(len_blpr, len_delasso, len_ri_M))
  # print(Cover)
  # print(">> Coverage")
  # print(colMeans(Cover))
  # print(">> Length")
  # print(colMeans(Len))
# }



nonIntersectPerm <- function(p){
    ind <- sample(p)
    while(any(ind == 1:p)){
      ind <- sample(p)
    }
    return(diag(length(ind))[ind,])
}

a = nonIntersectPerm(4)
print(a)





