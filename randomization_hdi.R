pacman::p_load(CVXR, fastclime, glmnet, hdi, mvtnorm, RPtests)


split_perm <- function(n, n_g, X){
  G <- list()
  XGX <- list()
  fir <- sample(n, n/2)
  sec <- setdiff(1:n, fir)
  for (i in 1:n_g){
    ind <- rep(0, n)
    ind[fir] <- sample(sec)
    ind[sec] <- sample(fir)
    # ind <- c(sample(sec), sample(fir))
    G[[i]] <- diag(length(ind))[ind,]
    XGX[[i]] <- (t(X) %*% G[[i]] %*% X) / n
  }
  return(list('G'=G, 'XGX'=XGX))
}


equal_flips <- function(n, n_g, X){
  G <- list()
  XGX <- list()
  ind <- 1:n
  signs <- rep(0, n)
  for (i in 1:n_g){
    fir <- sample(n, n/2)
    sec <- setdiff(ind, fir)
    signs[fir] <- 1
    signs[sec] <- -1
    G[[i]] <- diag(signs)
    XGX[[i]] <- (t(X) %*% G[[i]] %*% X) / n
  }
  return(list('G'=G, 'XGX'=XGX))
}


edit.fastclime.selector <- function(lambdamtx, icovlist, lambda){
  # Edited fastclime.selector to 
  # 1) ensure that ||M_{, i}||_1 is larger than tol
  # 2) include feasibility as output
  
  gcinfo(FALSE)
  d<-dim(icovlist[[1]])[2]
  maxnlambda<-dim(lambdamtx)[1]
  icov<-matrix(0,d,d)
  adaj<-matrix(0,d,d)
  seq<-rep(0,d)
  threshold<-1e-5
  status<-0
  
  for(i in 1:d){
    temp_lambda<-which(lambdamtx[,i]>lambda)
    seq[i]<-length(temp_lambda)
    if((seq[i]+1)>maxnlambda){
      status <- 1
      icov[,i] < -icovlist[[seq[i]]][,i]
    }
    else{
      icov[,i]<-icovlist[[seq[i]+1]][,i]
    }
  }
  
  icov <- icov*(abs(icov)<= abs(t(icov)))+ t(icov)*(abs(icov)> abs(t(icov)))
  
  tmpicov<-icov
  diag(tmpicov)<-0
  adaj<-Matrix(tmpicov>threshold, sparse=TRUE)*1
  
  sparsity<-(sum(adaj))/(d^2-d)
  
  if(status==1)
  {
    cat("Some columns do not reach the required lambda!\n You may want to increase lambda.min or use a larger nlambda. \n")
  }
  
  rm(temp_lambda,seq,d,threshold)
  gc()
  
  result<-list("icov"=icov, "adaj"=adaj,"sparsity"=sparsity, 'status'=status)
  class(result)="fastclime.selector"
  
  return(result)
}


bounded.fastclime.selector <- function(lambdamtx, icovlist, lambda, tol){
  # Edited fastclime.selector to 
  # 1) ensure that ||M_{, i}||_1 is larger than tol
  # 2) include feasibility as output
  
  gcinfo(FALSE)
  d<-dim(icovlist[[1]])[2]
  maxnlambda<-dim(lambdamtx)[1]
  icov<-matrix(0,d,d)
  adaj<-matrix(0,d,d)
  seq<-rep(0,d)
  threshold<-1e-5
  status<-0
  
  for(i in 1:d){
    temp_lambda<-which(lambdamtx[,i]>lambda)
    seq[i]<-length(temp_lambda)
    index <- 0
    if((seq[i]+1)>maxnlambda){
      status <- 1
      while (sum(abs(icovlist[[seq[i]-index]][,i])) > tol){
        index <- index + 1  
      }
      icov[,i] < -icovlist[[seq[i]-index]][,i]
    }
    else{
      while (sum(abs(icovlist[[seq[i]+1-index]][,i])) > tol){
        index <- index + 1  
      }
      icov[,i]<-icovlist[[seq[i]+1-index]][,i]
    }
  }
  
  icov <- icov*(abs(icov)<= abs(t(icov)))+ t(icov)*(abs(icov)> abs(t(icov)))
  
  tmpicov<-icov
  diag(tmpicov)<-0
  adaj<-Matrix(tmpicov>threshold, sparse=TRUE)*1
  
  sparsity<-(sum(adaj))/(d^2-d)
  
  if(status==1)
  {
    cat("Some columns do not reach the required lambda!\n You may want to increase lambda.min or use a larger nlambda. \n")
  }
  
  rm(temp_lambda,seq,d,threshold)
  gc()
  
  result<-list("icov"=icov, "adaj"=adaj,"sparsity"=sparsity, 'status'=status)
  class(result)="fastclime.selector"
  
  return(result)
}


rr_clime = function(y, X, n_g, n_solve, a_ind, a_val, n_ind, n_val, g_design, scale_res) {
  # M = individual rows of M^T obtained via solving the CLIME problem
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
  for (i in 1:t) { A[i, test_ind[i]] <- 1 }
  
  # Solve for M
  # Constants
  # X_bar <- colMeans2(X)
  # S <- t(X - X_bar) %*% (X - X_bar) / n
  S <- t(X) %*% (X) / n
  lambda <- sqrt(log(p) / n)
  
  # # CVXR/Gurobi
  # # Setup
  # M_ <- Variable(p, p)
  # obj <- Minimize(mixed_norm(M_, 1, 1))
  # 
  # # Define constraints
  # constraints <- list()
  # # for (i in 1:t){
  # for (i in 1:p){
  #   i_t <- rep(0, p)
  #   # i_t[[test_ind[i]]] <- 1
  #   i_t[i] <- 1
  #   # constraints <- c(constraints, S %*% t(M_[test_ind[i], ]) - i_t <= lambda)
  #   # constraints <- c(constraints, S %*% t(M_[test_ind[i], ]) - i_t >= -lambda)
  #   constraints <- c(constraints, S %*% t(M_[i, ]) - i_t <= lambda)
  #   constraints <- c(constraints, S %*% t(M_[i, ]) - i_t >= -lambda)
  # }
  # 
  # # Solve
  # prob <- Problem(obj, constraints)
  # result <- solve(prob, solver = "GUROBI")
  # M <- result$getValue(M_)
  
  # paralp
  # Setup
  M <- matrix(0, nrow=p, ncol=p)
  paralp_A <- rbind(cbind(S, -S), cbind(-S, S))
  paralp_c <- rep(-1, 2 * p)
  paralp_cbar <- rep(0, 2 * p)
  paralp_bbar <- rep(1, 2 * p)

  # Solve
  # for (i in 1:t){
  for (i in 1:p){
    paralp_b <- rep(0, 2 * p)
    # paralp_b[test_ind[i]] <- 1
    # paralp_b[(p+test_ind[i])] <- -1
    paralp_b[i] <- 1
    paralp_b[(p+i)] <- -1
    M_ <- paralp(paralp_c, paralp_A, paralp_b, paralp_cbar, paralp_bbar, lambda)
    # M[, test_ind[i]] <-  t(M_[1:p] - M_[(p+1):(2*p)])
    M[, i] <-  t(M_[1:p] - M_[(p+1):(2*p)])
  }

  M <- t(M)

  # Sqrt LASSO + correction
  sqrt_lasso = RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps = y - X %*% sqrt_lasso$beta
  beta_dlasso = A %*% sqrt_lasso$beta + 1 / n * A %*% M %*% t(X) %*% eps
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps = eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }
  
  An = (1 / sqrt(n)) * A %*% M %*% t(X)   # (s + # inactive) x n
  Tobs = t(sqrt(n) * (beta_dlasso - test_val))
  
  Tvals <- matrix(0, nrow=0, ncol=t)       # R x (s + # inactive)
  
  fir <- sample(n, n/2)
  sec <- setdiff(1:n, fir)
  
  for (i in 1:n_g){
    if(g_design=="perm"){
      ind <- rep(0, n)
      ind[fir] <- sample(sec)
      ind[sec] <- sample(fir)
      G <- diag(length(ind))[ind,]
    } 
    else if(g_design=="sign"){
      G <- diag(sample(c(rep(1, n/2), rep(-1, n/2))))
    }
    Tvals = rbind(Tvals, t(An %*% G %*% eps))
  }
  
  get_q = sapply(1:t, function(j){
    quantile(Tvals[,j], c(.025, .975)) / sqrt(n)
  })
  
  ci_a <- cbind(beta_dlasso[1:ac] - get_q[2,1:ac], beta_dlasso[1:ac] - get_q[1,1:ac])
  ci_n <- cbind(beta_dlasso[(ac+1):t] - get_q[2,(ac+1):t], beta_dlasso[(ac+1):t] - get_q[1,(ac+1):t])
  
  returnList <- list('ci_a' = ci_a, 'ci_n' = ci_n, 
                     'norm1_beta' = norm(as.matrix(sqrt_lasso$beta), '1'), 'norm1_dbeta' = norm(as.matrix(beta_dlasso), '1'),
                     'norm0_beta' = sum(abs(sqrt_lasso$beta) > 0), 'norm1_eps' = norm(as.matrix(eps), '1'))
  return(returnList)
}


rr_clime_all = function(y, X, n_g, M, a_ind, a_val, n_ind, n_val, g_design, scale_res) {
  # M = sparse M^T obtained via solving CLIME
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
  
  # Sqrt LASSO + correction
  beta_dlasso <- list()
  norm1_dbeta <- list()
  sqrt_lasso <- RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps <- y - X %*% sqrt_lasso$beta
  for (i in names(M)){
    beta_dlasso[[i]] <- sqrt_lasso$beta + 1 / n * M[[i]] %*% t(X) %*% eps
    norm1_dbeta <- norm(as.matrix(beta_dlasso[[i]]), '1')
  }
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps <- eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }
  
  An <- list()
  Tvals <- list()
  for (i in names(M)){
    An[[i]] <- (1 / sqrt(n)) * A %*% M[[i]] %*% t(X)   # (s + # inactive) x n
    Tvals[[i]] <- matrix(0, nrow=0, ncol=t)            # R x (s + # inactive)
  }
  # Tobs <- t(sqrt(n) * (A %*% beta_dlasso - test_val))
  
  # R = 1
  # while(R < n_g){
  #   if(g_design=="perm"){
  #     G = get_perm(n)
  #   } else if(g_design=="sign"){
  #     G = diag(sample(c(rep(1, n/2), rep(-1, n/2))))
  #   }
  #   if(max(abs(Matrix(G - diag(rep(1,n)))))>0){
  #     Tvals = rbind(Tvals, t(An %*% G %*% eps))
  #     R = R +1
  #   }
  # }
  
  # Get Tvals
  fir <- sample(n, n/2)
  sec <- setdiff(1:n, fir)
  
  for (i in 1:n_g){
    if(g_design=="perm"){
      ind <- rep(0, n)
      ind[fir] <- sample(sec)
      ind[sec] <- sample(fir)
      G <- diag(length(ind))[ind,]
    } 
    else if(g_design=="sign"){
      fir <- sample(n, n/2)
      sec <- setdiff(1:n, fir)
      ind <- rep(0, n)
      ind[fir] <- 1
      ind[sec] <- -1
      G <- diag(ind)
    }
    for (i in names(M)){
      Tvals[[i]] <- rbind(Tvals[[i]], t(An[[i]] %*% G %*% eps)) 
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
