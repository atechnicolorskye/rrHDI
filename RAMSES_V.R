library(Rcpp)
library(RcppArmadillo)


## Panos Toulis
## Implements the randomization method for regression models.
## Randomization Method for Standard ErrorS  (RAMSES)
get_clustered_eps = function(y, X, lam, lam0, clustering) {
  bhat = OLS_c(y, X)
  Q = matrix(lam, ncol=1)
  bhat_r = restricted_OLS_c(y, X, bhat, Q=Q, c=lam0)
  er = y - X %*% bhat_r # restricted residuals.
  stopifnot(all(!is.na(er))) ## no NAs in residuals.
  # return
  lapply(clustering, function(i) er[i])
}

one_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  srt = sort(tvals)
  M = length(tvals)
  k = ceiling(M * (1-alpha))
  Tk = srt[k]
  if(abs(tobs - Tk) < tol) {
    # if tobs=Tk
    ax = (M * alpha - sum(tvals > Tk)) / sum(abs(tvals - Tk) < tol)
    return(runif(1) <= ax) ## randomize.
  }
  
  return(tobs > Tk)
}

two_sided_test = function(tobs, tvals, alpha=0.05) {
  m1 = one_sided_test(tobs, tvals, alpha=alpha/2)
  m2 = one_sided_test(-tobs, -tvals, alpha=alpha/2)
  return(m1+m2)
}

CHECK_two_sided = function(nreps=1e4) {
  # tvals = c(-1, 1, 1, 1.5, 2.5, 3, 3.5)
  tvals = rnorm(20)
  # y = sapply(tvals, function(i) two_sided_test(i, tvals, 0.05))
  t0 = proc.time()[3]
  y = replicate(nreps, {
    i = sample(tvals, size=1)
    two_sided_test(i, tvals, 0.05)
  })
  t1 = proc.time()[3]
  print(sprintf("Time elapsed =%.2f secs, mean rej%% = %.2f%%", t1-t0, mean(y)*100))
  
  t0 = proc.time()[3]
  y1 = replicate(nreps, {
    i = sample(tvals, size=1)
    
    out_pval(list(tobs=i, tvals=tvals), F)
    
  })
  t1 = proc.time()[3]
  print(sprintf("Time elapsed =%.2f secs, mean rej%% = %.2f%%", t1-t0, mean(y)*100))
  print(table(y, y1))
}


out_pval = function(rtest_out, ret_pval) {
  ## tranforms the output of rtest_c to a p-value. 
  ## out = list(tobs, tvals)
  #
  tobs = rtest_out$tobs
  tvals = c(rtest_out$tvals)

  n_all = length(tvals)
  n_higher = sum(tvals > (tobs + 1e-12))
  n_lower = sum(tvals < (tobs - 1e-12))
  n_equal = n_all - n_lower - n_higher
 
  p1 = (n_equal + n_higher) / n_all  # P(T >= Tobs)
  p2 = (n_equal + n_lower) / n_all  # P(T <= Tobs)
  
  pval = min(p1, p2)
  if(ret_pval) return(pval)
  
  # 
  # OLD: return(pval <= 0.025)
  return(two_sided_test(tobs, tvals))  # this test is less conservative.
}


# riTest_exact = function(model, num_R=999, ret_pval=F) {
#   # Exact test for SLR
#   y = model$y; x= model$x; beta0 = model$beta0; beta1= model$beta1
#   
#   out = r_test_exact(y, x, beta0, beta1, num_R)
#   out_pval(out, ret_pval)
# }
  
riTest_P = function(model, clustering, num_R=999, ret_pval=F) {
  ## Randomization Test with Permutation only.
  # clustering = list(v1, v2, ...vk)  where vi partition {1, 2...n}. Should be ordered.
  # e.g., clustering = list((1, 2, 3), (4, 5), (6, 7, 8), (9, 10))
  #
  # model = y, X, lam, lam0
  #
  y = model$y; X = model$X; lam = model$lam; lam0 = model$lam0
  
  cl_eps_r = get_clustered_eps(y, X, lam, lam0, clustering)
  out = r_test_c(y, X, lam=lam, lam0 = lam0, cluster_eps_r=cl_eps_r, use_perm = T, use_sign = F, num_R = num_R)
  out_pval(out, ret_pval)
}

riTest_S = function(model, clustering, num_R=999, ret_pval=F) {
  ## Randomization Test with Sign-switch only.
  #
  y = model$y; X = model$X; lam = model$lam; lam0 = model$lam0
  
  cl_eps_r = get_clustered_eps(y, X, lam, lam0, clustering)
  if(length(clustering)==1) warning("Sign tests should have #clusters > 1.")
    
  out = r_test_c(y, X, lam=lam, lam0 = lam0, cluster_eps_r=cl_eps_r, use_perm = F, use_sign = T, num_R = num_R)
  out_pval(out, ret_pval)
}

riTest_Double = function(model, clustering, num_R=999, ret_pval=F) {
  ## Randomization Test with Permutation + Sign Switch
  #
  y = model$y; X = model$X; lam = model$lam; lam0 = model$lam0
  
  cl_eps_r = get_clustered_eps(y, X, lam, lam0, clustering)
  if(length(clustering)==1) warning("Sign tests should have #clusters > 1.")
  
  out = r_test_c(y, X, lam=lam, lam0 = lam0, cluster_eps_r=cl_eps_r, use_perm = T, use_sign = T, num_R = num_R)
  out_pval(out, ret_pval)
}

riTest_G = function(model, g_inv, num_R=999, ret_pval=F, ret_testout=F) {
  # y = vector of outcomes, X=covariates (first column should be ones)
  # num_R = #resamplings for the randomization method.
  # 
  # g_inv(e) = function of invariance for residuals
  #
  y = model$y; X = model$X; lam = model$lam; lam0 = model$lam0
  if(!all(X[, 1]==1)) {
    warning("No intercept.")
  }
  stopifnot(length(y)==nrow(X))
  stopifnot(length(lam)==ncol(X))
  
  Tn = function(eps) {
    b = fastLmPure(X+0, eps)$coefficients
    sum(lam * b)  # tn(ε) = λ' (X^t X)^-1  X^t ε
  }
  
  # 1. tobs
  bhat = fastLmPure(X+0, y)$coefficients
  tobs = (sum(lam * bhat) - lam0) ## = (t(X) * X)^-1 X^t e
  # same as tobs = Tn(y) - lam0
  
  # 2. e_r : restricted residuals.
  Q = matrix(lam, ncol=1)
  bhat_r = restricted_OLS_c(y, X, bhat, Q=Q, c=lam0)
  er = y - X %*% bhat_r # restricted residuals.
  tvals = c()
  for(r in 1:num_R) {
    # 1. iid
    er_new = g_inv(er)
    tvals = c(tvals, Tn(er_new))
  }
  
  out = list(tobs=tobs, tvals=tvals)
  ## Return the entire test output?
  if(ret_testout) return(out)
  ## Return pvalue or decision
  out_pval(out, ret_pval)
}

example = function(nreps) {
  # Example use of r_test.
  # Mean rejection should converge to 0.05
  #
  n = 100

  rej = c();
  rej1 = c(0)
  t0 = proc.time()[3]
  for(j in 1:nreps) {
    x = rnorm(n, mean=2)
    y = 1.0 + 1.0 * x + abs(2*x)*rnorm(n, sd=1)
    X = cbind(rep(1, n), x)
    
    lam = c(0, 1)
    lam0 = 1.0 + 0.0
    cl = list(1:n) # s.list(er) # list(clust1) #, clust2, clust3)
    model = list(y=y, X=X, lam=lam, lam0=lam0)
    # Test beta1=1
    # dec = riTest_P(model, clustering = cl)
    # dec = riTest_G(model, function(e) { sample(e) })
    dec = riTest_G(model, function(e) {
      g=sample(c(-1, 1), replace = T, size = n)
      g*e
    })
    rej = c(rej, dec)
    t1 = proc.time()[3]
    if(t1 - t0 > 5) {
      print(sprintf("iter= %d/%d -- rej = %.4f -- prev = %.3f", j, nreps, mean(rej), mean(rej1)))
      t0 = t1
    }
  }
}

# 
# st <- proc.time()
# example(10000) # JLK: about 10 times faster
# # example_fast(500) # JLK: with 500 nreps, mean rejection gets close to 0.05 almost every time
# et <- proc.time()
# print(sprintf("Fastest example took %f sec", (et-st)[3]))


# st <- proc.time()
# example(100) # JLK: original example
# et <- proc.time()
# print(sprintf("Original example took %f sec", (et-st)[3]))
