pacman::p_load(glmnet, hdi, mvtnorm, RPtests)

get_perm <- function(p){
  ind <- sample(p)
  while(any(ind == 1:p)){
    ind <- sample(p)
  }
  return(diag(length(ind))[ind,])
}


rr_dantzig = function(y, X, n_g, M, a_ind, a_val, n_ind, n_val, g_design, col_norm, scale_res) {
  # M = sparse M^T obtained via Dantzig selector
  # Test H0_j: beta_j = val_j individually
  
  n = nrow(X)
  p = ncol(X)
  stopifnot(length(y)==nrow(X))
  
  ac = length(a_ind)
  na = length(n_ind)
  t = ac + na
  
  test_ind = c(a_ind, n_ind)
  test_val = c(a_val, n_val) / col_norm[test_ind]
  
  Y = matrix(y, nrow=n, ncol=t, byrow=F) # n x p
  # Construct binary indicator for H0_j
  A = matrix(0, nrow=t, ncol=p)          # s x p
  for(i in 1:t) { A[i, test_ind[i]] = 1 }
  
  # Sqrt LASSO + correction
  sqrt_lasso = RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps = y - X %*% sqrt_lasso$beta
  beta_dlasso = sqrt_lasso$beta + 1 / n * M %*% t(X) %*% eps
  # eps = eps / col_norm
  # # Potentially wrong
  # if (scale_res == 1){
  #   # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
  #   eps = eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  # }
  
  An = (1 / sqrt(n)) * A %*% M %*% t(X)   # (s + # inactive) x n
  Tobs = t(sqrt(n) * (A %*% beta_dlasso - test_val))

  Tvals = matrix(0, nrow=0, ncol=t)      # R x (s + # inactive)
  
  R = 1
  while(R < n_g){
    if(g_design=="perm"){
      G = get_perm(n)
    } else if(g_design=="sign"){
      G = diag(sample(c(rep(1, n/2), rep(-1, n/2))))
    }
    # Check for self maps
    if(max(abs(Matrix(G - diag(rep(1,n)))))>0){
      Tvals = rbind(Tvals, t(An %*% G %*% eps))
      R = R +1
    }
  }
  
  get_q = sapply(1:t, function(j){
    quantile(Tvals[,j], c(.025, .975)) / sqrt(n)
  })
  
  # ci_a = cbind((beta_dlasso[a_ind] - get_q[2,1:ac]) / col_norm[a_ind], (beta_dlasso[a_ind] - get_q[1,1:ac]) / col_norm[a_ind])
  # ci_n = cbind((beta_dlasso[n_ind] - get_q[2,(ac+1):t]) / col_norm[n_ind], (beta_dlasso[n_ind] - get_q[1,(ac+1):t]) / col_norm[n_ind])
  ci_a = cbind(beta_dlasso[a_ind] - get_q[2,1:ac], beta_dlasso[a_ind] - get_q[1,1:ac])
  ci_n = cbind(beta_dlasso[n_ind] - get_q[2,(ac+1):t], beta_dlasso[n_ind] - get_q[1,(ac+1):t])
  
  returnList <- list('ci_a' = ci_a, 'ci_n' = ci_n, 
                     'norm1_beta' = norm(as.matrix(sqrt_lasso$beta / col_norm), '1'), 'norm1_dbeta' = norm(as.matrix(beta_dlasso / col_norm), '1'),
                     'norm0_beta' = sum(abs(sqrt_lasso$beta / col_norm) > 0), 'norm1_eps' = norm(as.matrix(eps / col_norm), '1'))
  return(returnList)
}


sel_M <- function(S, precisions, p, tol=2e-3){
  mu = rep(Inf, length(precisions))
  # Finds smallest mu
  for (i in 1:length(precisions)){
    mu[i] <- norm(diag(p) - S %*% precisions[[i]], 'M')
  }
  mu_star = min(mu) + tol
  idx = which.min(mu)
  # Finds M with smallest 1-norm that is mu + tol away from smallest mu
  for (i in (idx-1):1){
    if (mu[i] > mu_star){
      idx = i+1
      break
    }
  }
  norm_M_star = norm(precisions[[idx]], '1')
  returnList <- list("mu" = mu[idx], "mu_star" = mu_star,
                     'norm_M_star' = norm_M_star, 'M_star' = precisions[[idx]])
  return(returnList)
}


rr_ridge = function(y, X, n_g, a_ind, a_val, n_ind, n_val, g_design, scale_res){
  # Test H0_j: beta_j = val_j individually

  n = nrow(X)
  p = ncol(X)
  stopifnot(length(y)==nrow(X))

  ac = length(a_ind)
  na = length(n_val)
  t = ac + na

  test_ind = c(a_ind, n_ind)
  test_val = c(a_val, n_val)

  Y = matrix(y, nrow=n, ncol=t, byrow=F) # n x p
  # Construct binary indicator for H0_j
  A = matrix(0, nrow=t, ncol=p)          # s x p
  for(i in 1:t) { A[i, test_ind[i]] = 1 }

  # Sqrt LASSO + epsilon inflation
  sqrt_lasso = RPtests::sqrt_lasso(X, c(y), output_all = TRUE)
  eps = (y - X %*% sqrt_lasso$beta)
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps = eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }

  # Ridge + correction
  out = solve_ridge(y, X)
  beta_dridge = out$beta + out$lambda * out$P_inv %*% sqrt_lasso$beta
  # An = 1 / sqrt(n) * A %*% out$P_inv %*% t(X)   # (s + # inactive) x n
  An = sqrt(n) * A %*% out$P_inv %*% t(X)   # (s + # inactive) x n
  
  # # Ridge + projection correction
  # lambda = 1 / (10 * n)
  # I = diag(p)
  # P_inv = solve((t(X) %*% X) / n + lambda * I)
  # beta_ridge = P_inv %*% (t(X) %*% y) / n
  # P_X = t(X) %*% solve(X %*% t(X)) %*% X
  # beta_dridge = beta_ridge - ((matrix(1, p, p) - I) * P_X) %*% sqrt_lasso$beta
  # An = 1 / sqrt(n) * A %*% P_inv %*% t(X)   # (s + # inactive) x n

  Tobs = t(sqrt(n) * (A %*% beta_dridge - test_val))

  Tvals = matrix(0, nrow=0, ncol=t)         # R x (s + # inactive)

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
    quantile(Tvals[,j], c(.025, .975)) / sqrt(n)
  })

  ci_a = cbind(beta_dridge[a_ind]-get_q[2,1:ac], beta_dridge[a_ind]-get_q[1,1:ac])
  ci_n = cbind(beta_dridge[n_ind]-get_q[2,(ac+1):t], beta_dridge[n_ind]-get_q[1,(ac+1):t])

  returnList <- list("ci_a" = ci_a,"ci_n" = ci_n,
                     'norm1_beta' = norm(as.matrix(sqrt_lasso$beta), '1'), 'norm1_dbeta' = norm(as.matrix(beta_dridge), '1'),
                     'norm0_beta' = sum(abs(sqrt_lasso$beta) > 0), 'norm1_eps' = norm(as.matrix(eps), '1'))
  return(returnList)
}


solve_ridge = function(y, x) {
  stopifnot(length(y)==nrow(x))

  p = ncol(x)
  n = nrow(x)

  I = diag(p)
  
  # # Scale lambda_0, P_inv and beta w.r.t to glmnet objective function
  # lambda_0 = sqrt(log(p) / n) / n   # Equivalent to rescaling x and y by sqrt(n)
  # lambdas = seq(3, 0.1, by=-.1) * lambda_0
  # out = cv.glmnet(x, y, lambda=lambdas, nfolds=n, alpha=0, intercept=F)
  # lambda = out$lambda.min
  # P_inv = solve((t(x) %*% x) / n + lambda * I)
  # beta = P_inv %*% (t(x) %*% y) / n
  # beta_ <- as.vector(coef(out, x=x, y=y, s=out$lambda.min, exact = TRUE))[-1]
  # norm(beta - beta_, 'm)
  
  # browser()
  
  # Slightly different empirically but should technically be equivalent
  # https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  # Seems like not specifying lambda_0, and standardising x negatively impacts performance
  sd_y <- sqrt(var(y) * (n - 1) / n)[1, 1]
  # mean_x <- colMeans(x)
  # sd_x <- sqrt(apply(x, 2, var) * (n-1) / n)
  # x_sc <- matrix(NA, nrow = n, ncol = p)
  # for(i in 1:p){
  #   x_sc[,i] <- (x[,i] - mean_x[i]) / sd_x[i]
  # }
  lambda_0 = sqrt(log(p) /n) / n
  # lambda_0 = sqrt(log(p) /n)
  lambdas = seq(3, 0.1, by=-.1) * lambda_0
  # # out = cv.glmnet(sqrt(n) * x, sqrt(n) * y, lambda=lambdas, nfolds=n, alpha=0, intercept=F)
  out = cv.glmnet(x, y, lambda=lambdas, nfolds=n, alpha=0, intercept=F, standardize=F, thres=1e-8)
  lambda = out$lambda.min * n / sd_y
  # P_inv = solve((t(x_sc) %*% x_sc) + lambda * I)
  # beta = P_inv %*% (t(x_sc) %*% y)
  P_inv = solve((t(x) %*% x) + lambda * I)
  beta = P_inv %*% (t(x) %*% y)
  # beta_ <- as.vector(coef(out, x=x, y=y, s=out$lambda.min, exact = TRUE))[-1]
  
  list(beta=beta, lambda=lambda, P_inv=P_inv)
}


# rr_ridge_old = function(y, X, n_g, a_ind, a_val, n_ind, n_val, g_design, scale_res){
#   # Test H0_j: beta_j = val_j individually
# 
#   n = nrow(X)
#   p = ncol(X)
#   stopifnot(length(y)==nrow(X))
# 
#   ac = length(a_ind)
#   na = length(n_val)
#   t = ac + na
# 
#   test_ind = c(a_ind, n_ind)
#   test_val = c(a_val, n_val)
# 
#   Y = matrix(y, nrow=n, ncol=t, byrow=F) # n x p
#   # Construct binary indicator for H0_j
#   A = matrix(0, nrow=t, ncol=p)          # s x p
#   for(i in 1:t) { A[i, test_ind[i]] = 1 }
# 
#   # glmnet LASSO + epsilon inflation
#   beta_lasso = matrix(coef(cv.glmnet(X, y, intercept=F))[-1])
#   eps = y - X %*% beta_lasso
#   if (scale_res == 1){
#     # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
#     eps = eps * sqrt(n / (n - sum(abs(beta_lasso) > 0)))
#   }
# 
#   # Ridge + correction
#   out = solve_ridge(y, X)
#   beta_dridge = out$bhat + out$lam * out$P_inv %*% beta_lasso
#   An = 1 / sqrt(n) * A %*% out$P_inv %*% t(X)   # (s + # inactive) x n
# 
#   # # Ridge + projection correction
#   # lambda = 1 / (10 * n)
#   # I = diag(p)
#   # P_inv = solve((t(X) %*% X) / n + lambda * I)
#   # beta_ridge = P_inv %*% (t(X) %*% y) / n
#   # P_X = t(X) %*% solve(X %*% t(X)) %*% X
#   # beta_dridge = beta_ridge - ((matrix(1, p, p) - I) * P_X) %*% sqrt_lasso$beta
#   # An = 1 / sqrt(n) * A %*% P_inv %*% t(X)   # (s + # inactive) x n
# 
#   Tobs = sqrt(n) * t(A %*% beta_dridge - test_val)
# 
#   Tvals = matrix(0, nrow=0, ncol=t)             # R x (s + # inactive)
# 
#   R = 1
#   while(R < n_g){
#     if(g_design=="perm"){
#       G = get_perm(n)
#     } else if(g_design=="sign"){
#       G = diag(sample(c(rep(1, n/2), rep(-1, n/2))))
#     }
#     if(max(abs(Matrix(G - diag(rep(1,n)))))>0){
#       Tvals = rbind(Tvals, t(An %*% G %*% eps))
#       R = R +1
#     }
#   }
# 
#   get_q = sapply(1:t, function(j){
#     quantile(Tvals[,j], c(.025, .975)) / sqrt(n)
#   })
# 
#   ci_a = cbind(beta_dridge[a_ind]-get_q[2,1:ac], beta_dridge[a_ind]-get_q[1,1:ac])
#   ci_n = cbind(beta_dridge[n_ind]-get_q[2,(ac+1):t], beta_dridge[n_ind]-get_q[1,(ac+1):t])
# 
#   returnList <- list("ci_a" = ci_a,"ci_n" = ci_n,
#                      'norm1_beta' = norm(as.matrix(beta_lasso), '1'), 'norm1_dbeta' = norm(as.matrix(beta_dridge), '1'),
#                      'norm0_beta' = sum(abs(beta_lasso) > 0), 'norm1_eps' = norm(as.matrix(eps), '1'))
#   return(returnList)
# }
# 
# 
# solve_lasso= function(y, x) {
#   stopifnot(length(y) == nrow(x))
#   n = length(y)
#   x1 = scale(x,TRUE,TRUE)/sqrt(n-1)
#   #
#   out = lars(x1, y, intercept=F)
#   E  = apply(coef(out), 1, function(row) which(is.active(row)))
#   out1 = cv.glmnet(x1, y, intercept=F)
#   c1 = which(is.active(coef(out1))) -1
#   sel = which(unlist(lapply(E, function(l) setequal(l, c1))))
# 
#   if(length(sel) == 0) { return(list()) }
# 
#   # sel = selected coefficients based on CV.
#   # out = cv.lars(x, y)
#   bhat = coef(out)[sel,];
#   lam1 = out$lambda[sel]
#   # return the best-CV LASSO solution.
#   list(lambda=lam1, betahat=bhat, x1=x1, S=which(is.active(bhat)))
# }
# 
# 
# solve_ridge = function(y, x) {
#   stopifnot(length(y)==nrow(x))
#   out  = cv.glmnet(x, y, alpha=0, intercept=F)
#   lam = out$lambda.min
#   p = ncol(x)
#   I = diag(p)
#   P = t(x) %*% x + lam * I; P_inv = solve(P)
#   #
#   bhat = P_inv %*% t(x) %*% y
#   list(bhat=bhat, lam=lam, P=P, P_inv=P_inv)
# }
