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
nsim = 1
# b = 25
Cover = matrix(0,  nrow=0, ncol=3)
Len = matrix(0, nrow=0, ncol=3)

colnames(Cover) = c("BLPR", "DeLASSO", "Randomiz.")
colnames(Len) = c("BLPR", "DeLASSO", "Randomiz.")

# intlen = c()

for (i in 1:nsim) {
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
  
  ## residual bootstrap Lasso+Partial Ridge
  print("> Res. bootstrap + Ridge")
  out_lpr = bootLPR(x = x, y = y, type.boot = "residual", B = 500,
                    parallel.boot = TRUE, ncores.boot = 2)

  ci = t(out_lpr$interval.LPR)[ind,]
  # print(mean(ci[,2] - ci[,1]))
  stopifnot(nrow(ci)==length(beta_0))
  intv = (ci[,1] <= beta_0) & (ci[,2] >= beta_0)
  cover_blpr = mean(intv)  # mean coverage
  len_blpr = mean(ci[,2] - ci[,1])  # mean length.
  
  # 2. Debiased Lasso
  print("> Debiased LASSO")
  # out_lasso = SSLasso(x, y, intercept = FALSE)
  # ci = cbind(out_lasso$low.lim[ind], out_lasso$up.lim[ind])
  # stopifnot(nrow(ci)==length(beta_0))
  # intv = (ci[,1] <= beta_0) & (ci[,2] >= beta_0)
  # cover_delasso = mean(intv)
  # len_delasso = mean(ci[, 2] - ci[, 1])  # mean length.
  
  ## Randomization
  print("> Residual Randomization")
  
  # DantzigM
  S <- (t(x) %*% x) / n
  lambda_min <- min(0.1, 2 * sqrt(log(p) / n))
  n_solve <- 1000
  dantzig_M <- fastclime(S, lambda_min, n_solve)
  M <- fastclime.selector(dantzig_M$lambdamtx, dantzig_M$icovlist, 2 * sqrt(log(p) / n))
  out3 = riSingle_hdi_dantzig(y, x, M$icov, S_test=ind, lam0_vals=beta_0)
  cover_ri = mean(out3$pvals >=  0.025)
  len_ri = mean(out3$ci[, 2] - out3$ci[, 1]) # mean length.
  
  # Cover = rbind(Cover, c(cover_blpr, cover_delasso, cover_ri))
  # Len = rbind(Len, c(len_blpr, len_delasso, len_ri))
  # print(Cover)
  # print(">> Coverage")
  # print(colMeans(Cover))
  # print(">> Length")
  # print(colMeans(Len))
}

covprob = colSums(coverage)/nsim
avgintlen = colSums(intlen)/nsim

Results = as.data.frame(cbind(covprob, avgintlen))

# save(Results, file=sprintf("out/sim_blpr_sims%d_b%d.rda", nsim, b)

