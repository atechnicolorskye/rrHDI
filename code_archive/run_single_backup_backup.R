## RRI simulations for NeurIPS 2020 
rm(list=ls())
library("doParallel")
library("fastclime")
library("glmnet")
# library("lcmix")
library("mvtnorm")

library("HDCI") # Buhlmann
source("lasso_inference.r") # Montanari
source("randomization_hdi.R") # RRI

main_sim = function(s0, X_design, beta_design, err_design, g_design, nsim=100, n=50, p=200){
  
  print(sprintf("--===-  s0:          %d ==--===", s0))
  print(sprintf("--===-  X_design:    %s ==--===", X_design))
  print(sprintf("--===-  beta_design: %s ==--===", beta_design))
  print(sprintf("--===-  err_design:  %s ==--===", err_design))
  
  cols = c( "cov_blpr", "cov_dlasso", "cov_rri", 
            "len_blpr", "len_dlasso", "len_rri", 
            "s0", "Xs", "betas", "errs")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  Results = as.data.frame(Results)
  
  ### Now is the main Iteration.
  Q = matrix(0, nrow=0, ncol=6)  # interim results.
  colnames(Q) = colnames(Results)[1:6]
  
  rho = 0.8
  Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
  
  for (i in 1:nsim){
    print(sprintf("--===-  %d/%d iter ==--===",i, nsim))
    
    ######################### START DGP ############################
    ## Generate X
    x = NA
    if(X_design=="N1"){
      x = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma=diag(p))
    } else if(X_design=="G1"){
      x = rmvgamma(n, shape=1, rate=1, corr=diag(p))-1
    } else if(X_design=="N2"){
      x = matrix((sample(c(-2, 2), size = n*p, replace = T) + rnorm(n*p)), nrow = n, byrow = TRUE)
    } else if(X_design=="TG"){
      x = mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Tptz)
    } else if(X_design=="TGM"){
      x = rmvgamma(n, shape=1, rate=1, corr=Tptz)-1
    }
    
    ## Generate beta
    beta = rep(0, p)
    if(beta_design=="D1"){
      beta[1:s0] = sample(c(-1, 1), size=s0, replace=T)
    } else if(beta_design=="D5") {
      beta[1:s0] = sample(c(-5, 5), size=s0, replace=T)
    }
    # shuffle.
    beta = sample(beta)
    # focus on 15 largest
    ind = rev(order(abs(beta)))[1:15] # focus on those tests.
    beta_0 = beta[ind]
    
    ## Generate ERROR
    err = NA
    if(err_design=="N1"){
      err = rnorm(n)
    } else if (err_design=="G1"){
      err = rgamma(n, shape=1, rate = 1)
    } else if(err_design=="N2"){
      err = sample(c(-2, 2), size = n, replace = T) + rnorm(n)
    } else if(err_design=="HG"){
      err = mvtnorm::rmvnorm(n = 1, mean=rep(0, n), sigma = diag(rowSums(x*x)/p))
    } else if(err_design=="HMG"){
      err = sample(c(-2, 2), size = n, replace = T) + 
        mvtnorm::rmvnorm(n = 1, mean=rep(0, n), sigma = diag(rowSums(x*x)/p))
    }
    
    y = x %*% beta + as.vector(err)
    ######################### END DGP ############################
    
    ##################### START EXPERIMENT #######################

    ## 1. residual bootstrap Lasso+Partial Ridge
    print("> Res. bootstrap + Ridge")
    out_lpr = bootLPR(x = x, y = y, type.boot = "residual", B = 500, 
                      parallel.boot = TRUE, ncores.boot = 2)
    
    ci = t(out_lpr$interval.LPR)[ind,]
    stopifnot(nrow(ci)==length(beta_0))
    intv = (ci[,1] <= beta_0) & (ci[,2] >= beta_0)
    cover_blpr = mean(intv)  # mean coverage
    len_blpr = mean(ci[,2] - ci[,1])  # mean length.
    
    # 2. Debiased Lasso
    print("> Debiased LASSO")
    out_lasso = SSLasso(x, y, intercept = FALSE)
    ci = cbind(out_lasso$low.lim[ind], out_lasso$up.lim[ind])
    stopifnot(nrow(ci)==length(beta_0))
    intv = (ci[,1] <= beta_0) & (ci[,2] >= beta_0)
    cover_delasso = mean(intv)
    len_delasso = mean(ci[, 2] - ci[, 1])  # mean length.
    
    ## 3. Residual Randomization
    print("> Residual Randomization")
    # DantzigM
    S <- (t(x) %*% x) / n
    lambda_min <- min(.1, sqrt(log(p) / n))
    n_solve <- 500
    dantzig_M <- fastclime(S, lambda.min=lambda_min, nlambda=n_solve)
    M = dantzig_M$icovlist[[length(dantzig_M$icovlist)]]
    out3 = riSingle_hdi_dantzig(y, x, t(M), S_test=ind, lam0_vals=beta_0, g_design)
    
    # S <- (t(x) %*% x) / n
    # lambda_min <- min(0.1, 2 * sqrt(log(p) / n))
    # n_solve <- 1000
    # dantzig_M <- fastclime(S, lambda_min, n_solve)
    # M <- fastclime.selector(dantzig_M$lambdamtx, dantzig_M$icovlist, 2 * sqrt(log(p) / n))
    # out3 = riSingle_hdi_dantzig(y, x, t(M$icov), S_test=ind, lam0_vals=beta_0, g_design)
    
    cover_ri = mean(out3$pvals >=  0.025)
    len_ri = mean(out3$ci[, 2] - out3$ci[, 1]) # mean length.
    
    ############### Results ############### 
    Q = rbind(Q, c(cover_blpr, cover_delasso, cover_ri, len_blpr, len_delasso, len_ri))
    print(Q)
    print(">> Coverage and Length")
    print(colMeans(Q))
  }
  
  ## Replications over. Save results.
  print("--===-    RESULTS ==--===")
  cm = as.data.frame(t(colMeans(Q, na.rm=T)))
  R1 = cbind(cm, data.frame(s0=s0, Xs=X_design, betas=beta_design, errs=err_design, gs=g_design))
  Results = rbind(Results, R1)
  save(Results, file=sprintf("out/nsim%d_s%d_n%d_p%d_%s_%s_%s_%s.rda", 
                             nsim, s0, n, p, X_design, beta_design, err_design, g_design))
  print(Results)
  print("--===- --=-=== ==--===")

}

main_sim(10, X_design="N1", beta_design="D1", err_design="N1", g_design="perm", nsim=3, n=50, p=200)

# main_sim(5, X_design="N1", beta_design="D1", err_design="N1", g_design="sign", nsim=3, n=50, p=200)

# main_sim(5, X_design="TG", beta_design="D1", err_design="HG", g_design="sign", nsim=3, n=50, p=200)








