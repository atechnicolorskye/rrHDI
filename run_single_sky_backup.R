## RRI simulations for NeurIPS 2020 
rm(list=ls())
library(optparse)

option_list = list(
  make_option(c("-d", "--sparsity"), type="integer"),
  make_option(c("-x", "--x_design"), type="character"),
  make_option(c("-b", "--b_design"), type="character"),
  make_option(c("-e", "--e_design"), type="character"),
  make_option(c("-g", "--p_design"), type="character"),
  make_option(c("-r", "--nsim"), type="integer"),
  make_option(c("-n", "--n"), type="integer"),
  make_option(c("-p", "--p"), type="integer"),
  make_option(c("-s", "--nsolve"), type="integer")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opt = parse_args(opt_parser)

main_sim = function(s0, X_design, beta_design, err_design, g_design, nsim, n, p, nsolve){
  print(sprintf("--===-  s0:          %d ==--===", s0))
  print(sprintf("--===-  X_design:    %s ==--===", X_design))
  print(sprintf("--===-  beta_design: %s ==--===", beta_design))
  print(sprintf("--===-  err_design:  %s ==--===", err_design))
  
  cols = c( "cov_blpr", "cov_hdi", "cov_dlasso", "cov_rr", 
            "len_blpr", "len_hdi", "len_dlasso", "len_rr", 
            "s0", "Xs", "betas", "errs")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  Results = as.data.frame(Results)
  
  ### Now is the main Iteration.
  Q = matrix(0, nrow=0, ncol=8)  # interim results.
  colnames(Q) = colnames(Results)[1:8]
  
  rho = 0.8
  Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
  
  library(fastclime)
  library(glmnet)
  library(lcmix)
  library(mvtnorm)

  library(doParallel)

  library(HDCI) # Liu et. al
  library(hdi) # Buhlmann
  source("lasso_inference.R") # Montanari
  source("randomization_hdi.R") # RRI

  set.seed(2020)
  # registerDoParallel(4)

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
    out_blpr = bootLPR(x = x, y = y, type.boot = "residual", B = 500)
                      # parallel.boot = TRUE, ncores.boot = 2)

    ci_blpr = t(out_blpr$interval.LPR)[ind,]
    stopifnot(nrow(ci_blpr)==length(beta_0))
    intv_blpr = (ci_blpr[,1] <= beta_0) & (ci_blpr[,2] >= beta_0)
    cover_blpr = mean(intv_blpr)  # mean coverage
    len_blpr = mean(ci_blpr[,2] - ci_blpr[,1])

    ## 2. Buhlmann
    print("> Bulhmann")
    out_hdi = boot.lasso.proj(x, y, family = "gaussian", standardize = TRUE,
                              multiplecorr.method = "WY",
                              parallel = TRUE, ncores = getOption("mc.cores", 2L),
                              betainit = "cv lasso", sigma = NULL, Z = NULL, verbose = FALSE,
                              return.Z = FALSE, robust= FALSE,
                              B = 500, boot.shortcut = FALSE,
                              return.bootdist = TRUE, wild = TRUE,
                              gaussian.stub = FALSE)

    ci_hdi = confint(out_hdi, level = 0.95, parm=ind)
    stopifnot(nrow(ci_hdi)==length(beta_0))
    intv_hdi = (ci_hdi[,1] <= beta_0) & (ci_hdi[,2] >= beta_0)
    cover_hdi = mean(intv_hdi)
    len_hdi = mean(ci_hdi[,2] - ci_hdi[,1])

    # 2. Debiased Lasso
    print("> Debiased LASSO")
    out_dlasso = SSLasso(x, y, intercept = FALSE)

    ci_dlasso = cbind(out_dlasso$low.lim[ind], out_dlasso$up.lim[ind])
    stopifnot(nrow(ci_dlasso)==length(beta_0))
    intv_dlasso = (ci_dlasso[,1] <= beta_0) & (ci_dlasso[,2] >= beta_0)
    cover_dlasso = mean(intv_dlasso)
    len_dlasso = mean(ci_dlasso[,2] - ci_dlasso[,1])
    
    ## 3. Residual Randomization
    print("> Residual Randomization")
    # DantzigM
    S <- (t(x) %*% x) / n
    lambda_min <- min(.1, 2 * sqrt(log(p) / n))
    dantzig_M <- fastclime(S, lambda.min=lambda_min, nlambda=nsolve)
    M = dantzig_M$icovlist[[length(dantzig_M$icovlist)]]
    out_rri = riSingle_hdi_dantzig(y, x, t(M), S_test=ind, lam0_vals=beta_0, g_design)
    
    ci_rr = out_rri$ci
    stopifnot(nrow(ci_rr)==length(beta_0))
    intv_rr = (ci_rr[,1] <= beta_0) & (ci_rr[,2] >= beta_0)
    cover_rr = mean(intv_rr)
    len_rr = mean(ci_rr[,2] - ci_rr[,1])
    
    ############### Results ############### 
    Q = rbind(Q, c(cover_blpr, cover_hdi, cover_dlasso, cover_rr, len_blpr, len_hdi, len_dlasso, len_rr))
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

# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$p_design, opt$nsim, opt$n, opt$p, opt$nsolve)
main_sim(10, X_design="N1", beta_design="D1", err_design="N1", g_design="perm", nsim=1, n=50, p=100, nsolve=100)
# main_sim(10, X_design="N1", beta_design="D1", err_design="N1", g_design="perm", n_solve=300, nsim=10, n=100, p=300)
# main_sim(10, X_design="N1", beta_design="D1", err_design="N1", g_design="perm", n_solve=300, nsim=10, n=100, p=400)

# main_sim(10, X_design="TG", beta_design="D1", err_design="G1", g_design="perm", n_solve=300, nsim=10, n=100, p=200)
# main_sim(10, X_design="TG", beta_design="D1", err_design="G1", g_design="perm", n_solve=300, nsim=10, n=100, p=300)
# main_sim(10, X_design="TG", beta_design="D1", err_design="G1", g_design="perm", n_solve=300, nsim=10, n=100, p=400)

# main_sim(10, X_design="TG", beta_design="D1", err_design="HG", g_design="sign", n_solve=300, nsim=50, n=100, p=200)
# main_sim(10, X_design="TG", beta_design="D1", err_design="HG", g_design="sign", n_solve=300, nsim=50, n=100, p=300)
# main_sim(10, X_design="TG", beta_design="D1", err_design="HG", g_design="sign", n_solve=300, nsim=50, n=100, p=400)








