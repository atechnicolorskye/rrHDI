pacman::p_load(glmnet, hdi, mvtnorm, RPtests)


split_perm <- function(n, n_g, X){
  XG <- array(0, dim=c(dim(X)[2], dim(X)[1], n_g))
  XGX_n <- array(0, dim=c(dim(X)[2], dim(X)[2], n_g))
  
  fir <- sample(n, n/2)
  sec <- setdiff(1:n, fir)
  for (i in 1:n_g){
    ind <- rep(0, n)
    ind[fir] <- sample(sec)
    ind[sec] <- sample(fir)
    G <- diag(length(ind))[ind,]
    XG[, , i] <- t(X) %*% G
    XGX_n[, , i] <- norm((XG[, , i] %*% X) / n, 'I')
  }
  return(list('XG'=XG, 'mean_XGX_n_I'=mean(XGX_n)))
}


equal_flip <- function(n, n_g, X){
  XG <- array(0, dim=c(dim(X)[2], dim(X)[1], n_g))
  XGX_n <- array(0, dim=c(dim(X)[2], dim(X)[2], n_g))
  
  ind <- 1:n
  signs <- rep(0, n)
  for (i in 1:n_g){
    fir <- sample(n, n/2)
    sec <- setdiff(ind, fir)
    signs[fir] <- 1
    signs[sec] <- -1
    G <- diag(signs)
    XG[, , i] <- t(X) %*% G
    XGX_n[, , i] <- norm((XG[, , i] %*% X) / n, 'I')
  }
  return(list('XG'=XG, 'mean_XGX_n_I'=mean(XGX_n)))
}


min.fastclime.selector <- function(lambdamtx, icovlist, lambda, test_ind, mean_XGX_n_I, C){
  # Edited to minimise lambda_i + ||M_{, i}||_1 E_Q [|X^T G X / n|_/infty]
  gcinfo(FALSE)
  d <- dim(icovlist[[1]])[2]
  maxnlambda <- dim(lambdamtx)[1]
  icov <- matrix(0,d,d)
  adaj <- matrix(0,d,d)
  seq <- rep(0,d)
  threshold<-1e-5
  
  icovlist_a <- array(0, c(maxnlambda, d, d))
  for (i in 1:maxnlambda){
    icovlist_a[i, , ] <- icovlist[[i]]
  }
  
  for(i in 1:d){
    if (i %in% test_ind){
      seq[i] <- length(which(lambdamtx[ ,i] > lambda))
      if((seq[i]+1) > maxnlambda){
        i_lambdas <- lambdamtx[, i][1:seq[i]]
        err <- C * i_lambdas + rowMaxs(abs(icovlist_a[1:seq[i], , i])) * mean_XGX_n_I
      }
      else{
        i_lambdas <- lambdamtx[, i][1:(seq[i]+1)]
        err <- C * i_lambdas + rowMaxs(abs(icovlist_a[1:(seq[i]+1), , i])) * mean_XGX_n_I
      }
      icov[, i] <- icovlist[[which.min(err)]][, i]
    }
    else{
      temp_lambda <- which(lambdamtx[, i] > lambda)
      seq[i] <- length(temp_lambda)
      if((seq[i]+1) > maxnlambda){
        icov[, i] <- icovlist[[seq[i]]][, i]
      }
      else{
        icov[, i] <- icovlist[[(seq[i]+1)]][, i]
      }
    }
  }
  
  icov <- icov*(abs(icov)<= abs(t(icov)))+ t(icov)*(abs(icov)> abs(t(icov)))
  
  tmpicov<-icov
  diag(tmpicov)<-0
  adaj<-Matrix(tmpicov>threshold, sparse=TRUE)*1
  
  sparsity<-(sum(adaj))/(d^2-d)
  
  rm(temp_lambda,seq,d,threshold)
  gc()
  
  result<-list("icov"=icov, "adaj"=adaj,"sparsity"=sparsity)
  class(result)="fastclime.selector"
  
  return(result)
}


rr_min_clime = function(y, X, n_g, clime_M, lambda, a_ind, a_val, n_ind, n_val, g_design, scale_res) {
  # M = sparse M^T obtained via solving modified CLIME problem
  # Test H0_j: beta_j = val_j individually
  
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(y)==nrow(X))
  
  ac <- length(a_ind)
  na <- length(n_ind)
  t <- ac + na
  
  test_ind <- c(a_ind, n_ind)
  test_val <- c(a_val, n_val)
  
  Y <- matrix(y, nrow=n, ncol=t, byrow=F) # n x p
  
  # Construct binary indicator for H0_j
  A <- matrix(0, nrow=t, ncol=p)          # s x p
  for(i in 1:t) { A[i, test_ind[i]] <- 1 }
  
  if(g_design=="perm"){
    group_actions <- split_perm(n, n_g, X)
  }
  else if(g_design=="sign"){
    group_actions <- equal_flip(n, n_g, X)
  }
  
  M <- list()
  for (i in c('10000')){
    M[[i]] <- min.fastclime.selector(clime_M$lambdamtx, clime_M$icovlist, lambda, test_ind, group_actions$mean_XGX_n_I, strtoi(i))$icov
  }
  
  beta_dlasso <- list()
  norm1_dbeta <- list()
  # Sqrt LASSO + correction
  sqrt_lasso <- RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps <- y - X %*% sqrt_lasso$beta
  for (i in names(M)){
    beta_dlasso[[i]] <- sqrt_lasso$beta + 1 / n * M[[i]] %*% t(X) %*% eps
    norm1_dbeta[[i]] <- norm(as.matrix(beta_dlasso[[i]]), '1')
  }
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps <- eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }
  
  An <- list()
  Tvals <- list()
  for (i in names(M)){
    An[[i]] <- (1 / sqrt(n)) * A %*% M[[i]]   # (s + # inactive) x n
    Tvals[[i]] <- matrix(0, nrow=0, ncol=t)   # R x (s + # inactive)
  }
  
  for (i in 1:n_g){
    for (j in names(M)){
      Tvals[[j]] <- rbind(Tvals[[j]], t(An[[j]] %*% group_actions$XG[, , i] %*% eps))
    }
  }
  
  get_q <- list()
  for (i in names(M)){
    get_q[[i]] <- sapply(1:t, function(j){
      quantile(Tvals[[i]][,j], c(.025, .975)) / sqrt(n)
    })
  }
  
  ci_a <- list()
  ci_n <- list()
  for (i in names(M)){
    ci_a[[i]] <- cbind(beta_dlasso[[i]][a_ind] - get_q[[i]][2,1:ac], beta_dlasso[[i]][a_ind] - get_q[[i]][1,1:ac])
    ci_n[[i]] <- cbind(beta_dlasso[[i]][n_ind] - get_q[[i]][2,(ac+1):t], beta_dlasso[[i]][n_ind] - get_q[[i]][1,(ac+1):t])
  }
  
  returnList <- list('ci_a' = ci_a, 'ci_n' = ci_n, 
                     'norm1_beta' = norm(as.matrix(sqrt_lasso$beta), '1'), 'norm1_dbeta' = norm1_dbeta,
                     'norm0_beta' = sum(abs(sqrt_lasso$beta) > 0), 'norm1_eps' = norm(as.matrix(eps), '1'))
  return(returnList)
}


rr_hdi = function(y, X, n_g, a_ind, a_val, n_ind, n_val, g_design, scale_res) {
  # M = sparse M by solving for infinity norm
  # Test H0_j: beta_j = val_j individually
  
  n = nrow(X)
  p = ncol(X)
  stopifnot(length(y)==nrow(X))
  
  ac = length(a_ind)
  na = length(n_ind)
  t = ac + na
  
  test_ind = c(a_ind, n_ind)
  test_val = c(a_val, n_val)
  
  Y = matrix(y, nrow=n, ncol=t, byrow=F) # n x p
  # Construct binary indicator for H0_j
  A = matrix(0, nrow=t, ncol=p)          # s x p
  for(i in 1:t) { A[i, test_ind[i]] = 1 }

  # Create data structures
  S = t(X) %*% X / n
  
  # Sample group actions
  perms <- split_perm(n, n_g, X)
  G <- perms$G
  XGX <- perms$XGX
  
  lambda <- sqrt(log(p) / n) * matrix(1, p, 1)
  
  # Get M
  # Setup
  tau <- Variable(n_g)
  # tau <- Variable(1)
  M_ <- Variable(t, p)
  obj <- Minimize(mean(tau))
  # obj <- Minimize(tau)
  
  # Define constraints
  constraints <- list()
  for (i in 1:t){
    i_t <- rep(0, p)
    i_t[[test_ind[i]]] <- 1
    # constraints <- c(constraints, S %*% t(M_[i, ]) - i_t <= tau)
    # constraints <- c(constraints, S %*% t(M_[i, ]) - i_t >= -tau)
    for (j in 1:n_g){
      constraints <- c(constraints, (XGX[[j]] - S) %*% t(M_[i, ]) + i_t <= tau[j] * matrix(1, p, 1))
      constraints <- c(constraints, (XGX[[j]] - S) %*% t(M_[i, ]) + i_t >= -tau[j] * matrix(1, p, 1))
      # constraints <- c(constraints, XGX[[j]] %*% t(M_[i, ]) <= lambda * matrix(1, p, 1))
      # constraints <- c(constraints, XGX[[j]] %*% t(M_[i, ]) >= -lambda * matrix(1, p, 1))
    }
  } 
  
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "GUROBI", feastol = 1e-3, num_iter = 1e6)
  M <- result$getValue(M_)
  
  # Sqrt LASSO + correction
  sqrt_lasso = RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps = y - X %*% sqrt_lasso$beta
  beta_dlasso = A %*% sqrt_lasso$beta + 1 / n * M %*% t(X) %*% eps
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps = eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }

  An = (1 / sqrt(n)) * M %*% t(X)         # (s + # inactive) x n
  Tvals = matrix(0, nrow=0, ncol=t)       # R x (s + # inactive)
  for (i in 1:n_g){
    Tvals = rbind(Tvals, t(An %*% G[[i]] %*% eps))
  }
  
  Tobs = t(sqrt(n) * (beta_dlasso - test_val))
  
  get_q = sapply(1:t, function(j){
    quantile(Tvals[,j], c(.025, .975)) / sqrt(n)
  })
  
  ci_a = cbind(beta_dlasso[1:ac] - get_q[2,1:ac], beta_dlasso[1:ac] - get_q[1,1:ac])
  ci_n = cbind(beta_dlasso[(ac+1):t] - get_q[2,(ac+1):t], beta_dlasso[(ac+1):t] - get_q[1,(ac+1):t])
  
  returnList <- list('ci_a' = ci_a, 'ci_n' = ci_n, 
                     'norm1_beta' = norm(as.matrix(sqrt_lasso$beta), '1'), 'norm1_dbeta' = norm(as.matrix(beta_dlasso), '1'),
                     'norm0_beta' = sum(abs(sqrt_lasso$beta) > 0), 'norm1_eps' = norm(as.matrix(eps), '1'))
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
