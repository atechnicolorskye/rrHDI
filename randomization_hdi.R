pacman::p_load(glmnet, hdi, Rmosek, RPtests)


actionGenerator <- function(X, n_g0, n_g = n_g0, type = "perm"){
  ### Input:
  # X : data matrix
  # n_g0 : initial number of draws for G_0
  # n_g : number of draws to output for G
  # type : "perm" or "sign"
  
  ### Output:
  # XG : list of p x n  matrices which are X^T G
  # XGX_n : list of p x p matrices which are X^T G X
  # max_XGX_n_I: max infinity norm of X^T G Xs
  ### actionList : list of actions (either row permutations or signs)
  
  n <- nrow(X)
  
  XG <- list()
  XGX_n <- list()
  actionList <- list()
  XGX_n_I <- rep(0, n_g0)
  max_XGX_n_I <- 0 
  
  ind <- rep(0, n)
  
  for (i in 1:n_g0){
    
    if(type == "perm"){
      
      fir <- sample(n, n/2)
      sec <- setdiff(1:n, fir)
      ind[fir] <- sample(sec)
      ind[sec] <- sample(fir)
      
      XG[[i]] <-  t(X[order(ind), ])

    } else if(type == "sign") {
      
      fir <- sample(n, n/2)
      sec <- setdiff(1:n, fir)
      ind[fir] <- 1
      ind[sec] <- -1
      
      XG[[i]] <- t(X * ind)
      
    }
    
    XGX_n[[i]] <- XG[[i]] %*% X / n
    
    ## save the L_inf norm of X^TGX
    XGX_n_I[i] <- max(abs(XGX_n[[i]]))
    if (XGX_n_I[i] > max_XGX_n_I){
      max_XGX_n_I <- XGX_n_I[i]
    }
  }
  
  ## check if selection is made independent of M
  if (n_g < n_g0){
    ## pick the n_g smallest
    sorted <- order(XGX_n_I)
    selected <- sorted[1:n_g]
    XG <- XG[selected]
    XGX_n <- XGX_n[selected]
    mean_XGX_n_I <- mean(XGX_n_I[selected])
    # max_XGX_n_I <- XGX_n_I[selected[n_g]]
  } else {
    mean_XGX_n_I <- mean(XGX_n_I)
  }
  
  # return( list(XGX_n = XGX_n[selected], actionList = actionList[selected]) )
  return( list(XG = XG, XGX_n = XGX_n, mean_XGX_n_I = mean_XGX_n_I, max_XGX_n_I = max_XGX_n_I) )
}


min.fastclime.selector.opt <- function(lambdamtx, icovlist, lambda, test_ind, mean_XGX_n_I, delta){
  # Adapted from the fastclime package to minimize delta * lambda + ||M_{, i}||_1 E_Q [|X^T G X / n|_/infty]
  # Optimized for coverage

  ### Input:
  # lambdamtx : list of lambdas, output from fastclime
  # icovlist : list of precision matrices, output from fastclime
  # lambda: user-defined lambda
  # test_ind: indices to test
  # mean_XGX_n_I: mean infinity norm of X^T G X
  # delta: user-defined weighting parameter for lambda

  ### Output:
  # M : p x p matrix

  gcinfo(FALSE)
  d <- dim(icovlist[[1]])[2]
  maxnlambda <- dim(lambdamtx)[1]
  M <- matrix(0, d, d)
  seq <- rep(0, d)

  icovlist_a <- array(0, c(maxnlambda, d, d))
  for (i in 1:maxnlambda){
    icovlist_a[i, , ] <- icovlist[[i]]
  }

  for (i in test_ind){
    seq[i] <- length(which(lambdamtx[, i] > lambda))
    if((seq[i]+1) > maxnlambda){
      i_lambdas <- lambdamtx[, i][1:seq[i]]
      err <- delta * i_lambdas + rowMaxs(abs(icovlist_a[1:seq[i], , i])) * mean_XGX_n_I
    }
    else{
      i_lambdas <- lambdamtx[, i][1:(seq[i]+1)]
      err <- delta * i_lambdas + rowMaxs(abs(icovlist_a[1:(seq[i]+1), , i])) * mean_XGX_n_I
    }
    M[, i] <- icovlist[[which.min(err)]][, i]
  }

  rm(seq, d)
  gc()

  result<-list("M"=M)
  class(result)="fastclime.selector"

  return(result)
}


min.fastclime.selector <- function(lambdamtx, icovlist, lambda, test_ind, mean_XGX_n_I, delta){
  # Adapted from the fastclime package to minimize delta * lambda + ||M_{, i}||_1 E_Q [|X^T G X / n|_/infty]
  
  ### Input:
  # lambdamtx : list of lambdas, output from fastclime
  # icovlist : list of precision matrices, output from fastclime
  # lambda: user-defined lambda
  # test_ind: indices to test
  # mean_XGX_n_I: mean infinity norm of X^T G X
  # delta: user-defined weighting parameter for lambda
  
  ### Output:
  # M : p x p matrix
  
  gcinfo(FALSE)
  
  d <- dim(icovlist[[1]])[2]
  maxnlambda <- dim(lambdamtx)[1]
  M <- matrix(0, d, d)
  seq <- rep(0, d)

  icovlist_a <- array(0, c(maxnlambda, d, d))
  for (i in 1:maxnlambda){
    icovlist_a[i, , ] <- icovlist[[i]]
  }
  
  for(i in 1:d){
    if (i %in% test_ind){
      seq[i] <- length(which(lambdamtx[, i] > lambda))
      if((seq[i]+1) > maxnlambda){
        i_lambdas <- lambdamtx[, i][1:seq[i]]
        err <- delta * i_lambdas + rowMaxs(abs(icovlist_a[1:seq[i], , i])) * mean_XGX_n_I
      }
      else{
        i_lambdas <- lambdamtx[, i][1:(seq[i]+1)]
        err <- delta * i_lambdas + rowMaxs(abs(icovlist_a[1:(seq[i]+1), , i])) * mean_XGX_n_I
      }
      M[, i] <- icovlist[[which.min(err)]][, i]
    }
    else{
      temp_lambda <- which(lambdamtx[, i] > lambda)
      seq[i] <- length(temp_lambda)
      if((seq[i]+1) > maxnlambda){
        M[, i] <- icovlist[[seq[i]]][, i]
      }
      else{
        M[, i] <- icovlist[[(seq[i]+1)]][, i]
      }
    }
  }
  
  M <- M * (abs(M) <= abs(t(M))) + t(M) * (abs(M) > abs(t(M)))
  
  rm(temp_lambda, seq, d)
  gc()
  
  result<-list("M"=M)
  class(result)="fastclime.selector"
  
  return(result)
}


# min_ell_infinity <- function(XGX_n, S, g, m_i)
# {
#   d <- dim(XGX_n)[1]
#   p <- dim(XGX_n)[3]
#   # reshape XGX_n_m_S
#   XGX_n_ <- aperm(XGX_n, c(1, 3, 2))
#   XGX_n <- array(0, c(prod(dim(XGX_n)[1:2]), dim(XGX_n)[3]))
#   for (i in 1:d){
#     XGX_n[((i - 1) * p + 1):(i * p), ] <- XGX_n_[i, ,]
#   }
# 
#   rm(XGX_n_)
#   gc()
#   
#   # Create constraints for t
#   XGX_constraints <- array(0, c(prod(dim(XGX_n)[1]), d+1))
#   for (i in 1:d){
#     XGX_temp <- array(0, c(p, d+1))
#     XGX_temp[, i] <- 1
#     XGX_constraints[((i - 1) * p + 1):(i * p), ] <- XGX_temp
#   }
#   # S_constraints <- array(0, c(p, d+1))
#   # S_constraints[, d+1] <- 1
#   S_constraints <- array(0, c(p, 1))
#   S_constraints[, 1] <- 1
#   # Linear constraints in [x; t]
#   prob <- list(sense="min")
#   # prob$A <- rbind(cbind(XGX_n, -XGX_constraints),
#   #                 cbind(XGX_n, XGX_constraints),
#   #                 cbind(-S, -S_constraints),
#   #                 cbind(-S, S_constraints))
#   prob$A <- rbind(cbind(-S, -S_constraints),
#                   cbind(-S, S_constraints))
#   # Bound values for constraints
#   # prob$bc <- rbind(blc=c(array(rep(-Inf, dim(XGX_n)[1])), array(rep(0, dim(XGX_n)[1])), rep(-Inf, p), -g),
#   #                  buc=c(array(rep(0, dim(XGX_n)[1])), array(rep(Inf, dim(XGX_n)[1])), -g, rep(Inf, p)))
#   prob$bc <- rbind(blc=c(rep(-Inf, p), -g),
#                    buc=c(-g, rep(Inf, p)))
#   # Bound values for variables
#   # prob$bx <- rbind(blx=c(rep(-Inf, p), rep(0, d+1)),
#   #                  bux=rep(Inf, p+d+1))
#   prob$bx <- rbind(blx=c(rep(-Inf, p), rep(0, 1)),
#                    bux=rep(Inf, p+1))
#   # Coefficients for variables
#   # prob$c <- c(rep(0, p), 0.01 * rep(1 / d,  d), 0.99)
#   prob$c <- c(rep(0, p), 1)
#   # prob$sol <- list(bas=list())
#   # prob$sol$bas$xx <- c(m_i, rep(0, d+1))
# 
#   prob$iparam <- list(OPTIMIZER="OPTIMIZER_INTPNT")
#   # prob$iparam <- list(OPTIMIZER="OPTIMIZER_FREE_SIMPLEX")
#   r <- mosek(prob, list(verbose=1))
#   # return(list('m'=r$sol$itr$xx[1:p]))
#   m <- r$sol$itr$xx[1:p]
#   print(r$sol$itr$xx[p+1])
#   rm(prob, r)
#   gc()
#   
#   return(list('m'=m))
#   
# }

lp_clime <- function(S, base, lambda)
{
  ### Input:
  # S : empirical covariance matrix
  # base : standard basis for ith coordinate
  # lambda: constraint for |e_i - m S|_infty
  
  ### Output:
  # m : solution for ith coordinate
  # prob : RMosek problem
  # sol : RMosek solution
  
  p <- nrow(S)
  
  # construct first set of constraints
  S_constraints <- array(0, c(p, p+1))
  S_constraints[, p+1] <- 1
  
  # combine constraints in [x; t]
  prob <- list(sense="min")
  prob$A <- rbind(cbind(-S, -S_constraints),
                  cbind(-S, S_constraints),
                  cbind(diag(p), -diag(p), array(0, c(p, 1))),
                  cbind(diag(p), diag(p), array(0, c(p, 1)))
  )
  # bound values for constraints
  prob$bc <- rbind(blc=c(rep(-Inf, p), -base, rep(-Inf, p), rep(0, p)),
                   buc=c(-base, rep(Inf, p), rep(0, p), rep(Inf, p)))
  
  # bound values for variables
  prob$bx <- rbind(blx=c(rep(-Inf, p), rep(0, p),  lambda),
                   bux=c(rep(Inf, 2*p), lambda))
  
  # coefficients for variables
  prob$c <- c(rep(0, p), rep(1, p), 0)
  
  # use Simplex algorithm
  prob$iparam <- list(OPTIMIZER="OPTIMIZER_FREE_SIMPLEX")
  
  # solve
  res <- mosek(prob, list(verbose=1))
  
  return(list('m'=res$sol$bas$xx[1:p], 'prob'=prob, 'sol'=res$sol$bas))
}


lp_clime_mod <- function(prob, sol, lambda)
{
  ### Input:
  # prob : RMosek problem from output of lp_clime
  # sol : RMosek solution from output of lp_clime
  # lambda: constraint for |e_i - m S|_infty
  
  ### Output:
  # m : solution for ith coordinate
  # sol : RMosek solution of modified problem
  
  # use previous solution
  prob$sol <- list(bas=sol)
  
  # change lambda
  p2_plus_1 <- ncol(prob$bx)
  prob$bx[, p2_plus_1] <- lambda
  
  # solve
  res <- mosek(prob, list(verbose=1))
  
  return(list('m'=res$sol$bas$xx[1:(0.5 * (p2_plus_1  - 1))], 'sol'=res$sol$bas))
}


rr_min_clime = function(y, X, n_g, clime_M, lambda, test_ind, test_val, g_design, scale_res) {
  # Test H0_j: beta_j = val_j for a high-dimensional linear model i.e. y = X^T beta + epsilon
  
  ### Input:
  # y : observations
  # X : covariates
  # n_g : number of group actions
  # clime_M : solution of CLIME problem (if applicable)
  # lambda : constraint for |I - M S|_infty (if applicable)
  # test_ind : indices of beta to test
  # test_val : values of beta to test
  # g_design : choice of permutation or sign-flip matrices for inference
  # scale_res : to scale results using a heuristic (1) or not (0)
  
  ### Output:
  # ci : confidence intervals of test indices
  
  # get n and p 
  n <- nrow(X)
  p <- ncol(X)
  
  # get number of indices for testing
  # ac <- length(a_ind)
  # na <- length(n_ind)
  # t <- length(test_ind)
  
  # test_ind <- c(a_ind, n_ind)
  # test_val <- c(a_val, n_val)
  
  # sample group actions
  if(g_design=="perm"){
    group_actions <- actionGenerator(X, n_g, type = "perm")
  }
  else if(g_design=="sign"){
    group_actions <- actionGenerator(X, n_g, type = "sign")
  }
  
  # compute empirical covariance
  S <- t(X) %*% (X) / n
  
  # M <- min.fastclime.selector(clime_M$lambdamtx, clime_M$icovlist, lambda, test_ind, group_actions$mean_XGX_n_I, 10000)$M
  
  aux_cf <- function(multiplier, base_sol, lambda, group_actions=NULL, base_err=NULL){
    # Auxiliary error function determining quality of an m_i given lambda 
    # Assumes multiplier * lambda < 1 and ensures only feasible m_is would be accepted  
    
    ### Input:
    # multiplier: multiplier for lambda
    # base_sol: output from lp_clime
    # lambda: constraint for |e_i - m S|_infty
    # group_actions: sampled group actions
    # base_err: base error given lambda = 2 * lambda
    
    ### Output:
    # aux_err: auxiliary error given lambda = multiplier * lambda, accounts for
    #          base error if group_actions is NULL

    out <- lp_clime_mod(base_sol$prob, base_sol$sol, multiplier * lambda)
    # checks if sampled group actions are accounted for
    if (is.null(group_actions) == FALSE){
      # computes error for multiplier * lambda, ensures only m_is yielding less than base_err are accepted
      err <- multiplier * lambda + mean(sapply(group_actions$XGX_n, function(Z){max(abs(out$m %*% Z))})) # sum(abs(out$m)) * group_actions$mean_XGX_n_I 
      aux_err <- multiplier * lambda + (err > base_err) + (out$sol$prosta != "PRIMAL_AND_DUAL_FEASIBLE")
    }
    else{
      aux_err <- par * lambda + (out$sol$prosta != "PRIMAL_AND_DUAL_FEASIBLE")
    }
    return(aux_err)
  }
  
  M <- array(0, c(p, p))
  mXG <- array(0, c(length(test_ind), n_g, n))
  for (j in 1:length(test_ind)){
    I_i <- rep(0, p)
    I_i[test_ind[j]] <- 1
    lambda <- sqrt(log(p) / n)
    # base_sol <- lp_clime(S, I_i, lambda)
    base_sol <- lp_clime(S, I_i, 0.99)
    # base_err <- 2 * lambda + sum(abs(base_sol$m)) * group_actions$mean_XGX_n_I
    base_err <- 0.99 + sum(abs(base_sol$m)) * group_actions$mean_XGX_n_I
    # base_err <- sqrt(log(p) / n) + mean(apply(abs(array_prod(base_sol$m, group_actions$XGX_n)), MAR=1, max))
    # ptm <- proc.time()
    # pre-allocate for loop
    # coeff <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
    # err <- array(coeff * sqrt(log(p) / n), c(length(coeff)))
    # m_s <- array(0, c(length(coeff), p))
    # aux_err <- array(0, c(length(coeff)))
    # for (c in 1:length(coeff)){
    #   m_s[c, ] <- lp_tf_clime(p, base_sol$prob, base_sol$sol, coeff[c] * sqrt(log(p) / n))$m
    #   aux_err[c] <- coeff[c] * sqrt(log(p) / n) + mean(apply(abs(array_prod(m_s[c, ], group_actions$XGX_n)), MAR=1, max))
    # }
    # err[aux_err >= base_err] <- Inf
    par <- optimize(f = aux_cf, interval=c(1e-3, (0.99 / lambda)), base_sol = base_sol, lambda = lambda, group_actions = group_actions, base_err = base_err)$minimum
    # par <- optimize(f = aux_cf, interval=c(1e-3, 2), base_sol = base_sol, lambda = lambda, group_actions = group_actions, base_err = base_err)$minimum
    print(par)
    # M[, idx] <- m_s[which.min(err), ]
    # M[, test_ind[j]]
    M[, test_ind[j]] <- lp_clime_mod(base_sol$prob, base_sol$sol, par * lambda)$m
    # sorted <- order(sapply(group_actions$XGX_n, function(Z){max(abs(M[, test_ind[j]] %*% Z))}))
    # selected <- sorted[1:(n_g / 5)]
    # mXG[j, , ] <- t(sapply(group_actions$XG[c(selected)], function(Z){M[, test_ind[j]] %*% Z}))
    mXG[j, , ] <- t(sapply(group_actions$XG, function(Z){M[, test_ind[j]] %*% Z}))
    # print(proc.time() - ptm)
    rm(base_sol)
    gc()
  }

  M <- t(M)
  
  # Sqrt LASSO + correction
  sqrt_lasso <- RPtests::sqrt_lasso(X, c(y), output_all = TRUE, intercept = FALSE)
  eps <- y - X %*% sqrt_lasso$beta
  beta_dlasso <- sqrt_lasso$beta + 1 / n * M %*% t(X) %*% eps
  norm1_dbeta <- norm(as.matrix(beta_dlasso), '1')
  if (scale_res == 1){
    # Inflates epsilons according to proposed heuristic from eq. 24 of Zhang and Cheng
    eps <- eps * sqrt(n / (n - sum(abs(sqrt_lasso$beta) > 0)))
  }
  
  # M <- M[test_ind, ]
  # Tvals <- matrix(0, nrow=n_g, ncol=length(test_ind))   # R x (s + # inactive)
  # for (i in 1:n_g){
  #   Tvals[i, ] <- t((1 / sqrt(n)) * M %*% group_actions$XG[[i]] %*% eps)
  # }
    
  Tvals <- matrix(0, nrow=n_g, ncol=length(test_ind))   # R x (s + # inactive)
  for (j in 1:length(test_ind)){
    Tvals[, j] <- (1 / sqrt(n)) * mXG[j, , ] %*% eps
  }
  
  get_q <- sapply(1:length(test_ind), function(j){
    quantile(Tvals[, j], c(.025, .975)) / sqrt(n)
  })
  
  # ci_a <- cbind(beta_dlasso[a_ind] - get_q[2,1:ac], beta_dlasso[a_ind] - get_q[1,1:ac])
  # ci_n <- cbind(beta_dlasso[n_ind] - get_q[2,(ac+1):t], beta_dlasso[n_ind] - get_q[1,(ac+1):t])

  ci <- cbind(beta_dlasso[test_ind] - get_q[2, ], beta_dlasso[test_ind] - get_q[1, ])

  return(list('ci' = ci))
}
