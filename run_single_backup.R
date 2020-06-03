## RRI simulations for NeurIPS 2020 
rm(list=ls())
library(fastclime)
library(glmnet)
library(lcmix)
library(mvtnorm)
library(optparse)
library(parallel)

library(HDCI) # Liu et. al
library(hdi) # Buhlmann
source("lasso_inference.R") # Montanari
library(SILM) # Zhang and Guang
source("randomization_hdi.R") # RR

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


SimE.CI = function (X, Y, set, M = 1000, alpha = 0.95)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  Gram <- t(X) %*% X/n
  score.nodewiselasso = getFromNamespace("score.nodewiselasso",
                                         "hdi")
  if (p > floor(n/2)) {
    node <- score.nodewiselasso(X, wantTheta = TRUE, verbose = FALSE,
                                lambdaseq = "quantile", parallel = FALSE, ncores = 2,
                                oldschool = FALSE, lambdatuningfactor = 1)
    Theta <- node$out
  }
  else {
    Theta <- solve(Gram)
  }
  sreg <- scalreg(X, Y)
  beta.hat <- sreg$coefficients
  sigma.sq <- sum((Y - X %*% beta.hat)^2)/(n - sum(abs(beta.hat) > 0))
  # Ensure confidence intervals are not infinite
  if (is.infinite(sigma.sq)){
    sigma.sq <- sum((Y - X %*% beta.hat)^2)
  }
  beta.db <- beta.hat + Theta %*% t(X) %*% (Y - X %*% beta.hat)/n
  Omega <- diag(Theta %*% Gram %*% t(Theta)) * sigma.sq
  stat.boot.st <- stat.boot.nst <- rep(NA, M)
  for (i in 1:M) {
    e <- rnorm(n)
    xi.boot <- Theta[set, ] %*% t(X) %*% e * sqrt(sigma.sq)/sqrt(n)
    stat.boot.nst[i] <- max(abs(xi.boot))
  }
  crit.nst <- quantile(stat.boot.nst, alpha)
  up.nst <- beta.db[set] + crit.nst/sqrt(n)
  low.nst <- beta.db[set] - crit.nst/sqrt(n)
  band.nst <- rbind(low.nst, up.nst)
  result <- list(beta.db[set], band.nst)
  names(result) <- c("de-biased Lasso", "band.nst")
  return(result)
}

main_sim = function(s0, X_design, beta_design, err_design, g_design, nsim, n_list, p_list, nsolve){
  print(sprintf("--===-  s:           %d ==--===", s0))
  print(sprintf("--===-  X_design:    %s ==--===", X_design))
  print(sprintf("--===-  beta_design: %s ==--===", beta_design))
  print(sprintf("--===-  err_design:  %s ==--===", err_design))
  print(sprintf("--===-  g_design:    %s ==--===", g_design))
  print(sprintf("--===-  nsim:        %s ==--===", nsim))
  print(sprintf("--===-  n:           %s ==--===", n_list))
  print(sprintf("--===-  p:           %s ==--===", p_list))
  print(sprintf("--===-  nsolve:      %s ==--===", nsolve))
  
  for (i in 1:length(n_list)){
    n = n_list[i]
    p = p_list[i]
    
    if (nsolve < 0){
      # Name columns
      cols = c( "a_cov_rr", "a_len_rr", 
                "n_cov_rr", "n_len_rr", 
                "s", "x", "b", "e")
      # Cache results
      Q = matrix(0, nrow=0, ncol=4) 
      colnames(Q) = cols[1:4]
    }
    else{
      cols = c( "a_cov_blpr", "a_cov_hdi", "a_cov_dlasso", "a_cov_slim", "a_cov_rr",
                "a_len_blpr", "a_len_hdi", "a_len_dlasso", "a_len_slim", "a_len_rr",
                "n_cov_blpr", "n_cov_hdi", "n_cov_dlasso", "n_cov_slim", "n_cov_rr",
                "n_len_blpr", "n_len_hdi", "n_len_dlasso", "n_len_slim", "n_len_rr",
                "s", "x", "b", "e")
      Q = matrix(0, nrow=0, ncol=20)
      colnames(Q) = cols[1:20]
    }
    Results = matrix(0, nrow=0, ncol=length(cols))
    colnames(Results) = cols
    Results = as.data.frame(Results)
    
    rho = 0.8
    Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
    
    seed = 0
    n_cores = detectCores(all.tests = FALSE, logical = TRUE)
    
    for (i in 1:nsim){
      seed = seed + 1
      set.seed(seed)
      
      print(sprintf("--===-  %d/%d iter ==--===",i, nsim))
      ######################### START DGP ############################
      ## Generate X
      x = NA
      # print(n, p)
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
      # Shuffle.
      beta = sample(beta)
      ind = rev(order(abs(beta)))
      # Pick active variables 
      ind_0 =  ind[1:s0]
      beta_0 = beta[ind_0]
      # Pick 15 inactive variables
      ind_1 = ind[(s0+1):(s0+15)]
      beta_1 = beta[ind_1]
      
      ## Generate errors.
      err = NA
      if(err_design=="N1"){
        err = rnorm(n)
      } else if (err_design=="G1"){
        err = rgamma(n, shape=1, rate = 1)
      } else if(err_design=="N2"){
        err = sample(c(-2, 2), size = n, replace = T) + rnorm(n)
      } else if(err_design=="HG"){
        err = mvtnorm::rmvnorm(n = 1, mean=rep(0, n), sigma = 2 * diag(rowSums(x*x)/p))
      } else if(err_design=="HMG"){
        err = sample(c(-2, 2), size = n, replace = T) + 
          mvtnorm::rmvnorm(n = 1, mean=rep(0, n), sigma = 2 * diag(rowSums(x*x)/p))
      }
      
      y = x %*% beta + as.vector(err)
      ######################### END DGP ############################
      
      ##################### START EXPERIMENT #######################
      if (nsolve > 0) {
        ## 1. Residual Bootstrap Lasso + Partial Ridge, parallel somehow doesn't play nice
        print("> Res. Bootstrap + Ridge")
        out_blpr = bootLPR(x = x, y = y, type.boot = "residual", B = 1000,
                           parallel = TRUE, parallel.boot = TRUE, ncores.boot = n_cores)

        ci_blpr_a = t(out_blpr$interval.LPR)[ind_0,]
        ci_blpr_n = t(out_blpr$interval.LPR)[ind_1,]
        stopifnot(nrow(ci_blpr_a)==length(beta_0))
        stopifnot(nrow(ci_blpr_n)==length(beta_1))
        intv_blpr_a = (ci_blpr_a[,1] <= beta_0) & (ci_blpr_a[,2] >= beta_0)
        intv_blpr_n = (ci_blpr_n[,1] <= beta_1) & (ci_blpr_n[,2] >= beta_1)
        cov_blpr_a = mean(intv_blpr_a)
        cov_blpr_n = mean(intv_blpr_n)
        len_blpr_a = mean(ci_blpr_a[,2] - ci_blpr_a[,1])
        len_blpr_n = mean(ci_blpr_n[,2] - ci_blpr_n[,1])

        ## 2. Buhlmann
        print("> Bulhmann")
        # Use wild when g_design is sign
        if (g_design == "sign") {
          wild = TRUE
        }
        else{
          wild = FALSE
        }
        # print(wild)
        out_hdi = boot.lasso.proj(x, y, family = "gaussian", standardize = TRUE,
                                  multiplecorr.method = "WY",
                                  parallel = TRUE, ncores = n_cores,
                                  betainit = "cv lasso", sigma = NULL, Z = NULL, verbose = FALSE,
                                  return.Z = FALSE, robust= FALSE,
                                  B = 1000, boot.shortcut = FALSE,
                                  return.bootdist = TRUE, wild = wild,
                                  gaussian.stub = FALSE)

        ci_hdi_a = confint(out_hdi, level = 0.95, parm=ind_0)
        ci_hdi_n = confint(out_hdi, level = 0.95, parm=ind_1)
        stopifnot(nrow(ci_hdi_a)==length(beta_0))
        stopifnot(nrow(ci_hdi_n)==length(beta_1))
        # intv_hdi = (ci_hdi[,1] <= beta_0) & (ci_hdi[,2] >= beta_0)
        # cover_hdi = mean(intv_hdi)
        # len_hdi = mean(ci_hdi[,2] - ci_hdi[,1])
        intv_hdi_a = (ci_hdi_a[,1] <= beta_0) & (ci_hdi_a[,2] >= beta_0)
        intv_hdi_n = (ci_hdi_n[,1] <= beta_1) & (ci_hdi_n[,2] >= beta_1)
        cov_hdi_a = mean(intv_hdi_a)
        cov_hdi_n = mean(intv_hdi_n)
        len_hdi_a = mean(ci_hdi_a[,2] - ci_hdi_a[,1])
        len_hdi_n = mean(ci_hdi_n[,2] - ci_hdi_n[,1])

        ## 3. Debiased Lasso
        print("> Debiased LASSO")
        out_dlasso = SSLasso(x, y, intercept = FALSE)

        ci_dlasso_a = cbind(out_dlasso$low.lim[ind_0], out_dlasso$up.lim[ind_0])
        ci_dlasso_n = cbind(out_dlasso$low.lim[ind_1], out_dlasso$up.lim[ind_1])
        stopifnot(nrow(ci_dlasso_a)==length(beta_0))
        stopifnot(nrow(ci_dlasso_n)==length(beta_1))
        intv_dlasso_a = (ci_dlasso_a[,1] <= beta_0) & (ci_dlasso_a[,2] >= beta_0)
        intv_dlasso_n = (ci_dlasso_n[,1] <= beta_1) & (ci_dlasso_n[,2] >= beta_1)
        cov_dlasso_a = mean(intv_dlasso_a)
        cov_dlasso_n = mean(intv_dlasso_n)
        len_dlasso_a = mean(ci_dlasso_a[,2] - ci_dlasso_a[,1])
        len_dlasso_n = mean(ci_dlasso_n[,2] - ci_dlasso_n[,1])
        
        ## 4. SILM
        print("> SILM")
        ci_silm_a = matrix(0, s0, 2)
        ci_silm_n = matrix(0, 15, 2)
        for (i in 1:s0){
          ci_silm = SimE.CI(x, y, ind_0[i], M = 500, alpha = 0.95)$band.nst
          ci_silm_a[i, ] = t(ci_silm)
        }
        for (i in 1:15){
          ci_silm = SimE.CI(x, y, ind_1[i], M = 500, alpha = 0.95)$band.nst
          ci_silm_n[i, ] = t(ci_silm)
        }
        stopifnot(nrow(ci_silm_a)==length(beta_0))
        stopifnot(nrow(ci_silm_n)==length(beta_1))
        intv_silm_a = (ci_silm_a[,1] <= beta_0) & (ci_silm_a[,2] >= beta_0)
        intv_silm_n = (ci_silm_n[,1] <= beta_1) & (ci_silm_n[,2] >= beta_1)
        cov_silm_a = mean(intv_silm_a)
        cov_silm_n = mean(intv_silm_n)
        len_silm_a = mean(ci_silm_a[,2] - ci_silm_a[,1])
        len_silm_n = mean(ci_silm_n[,2] - ci_silm_n[,1])
      }
      
      # Residual Randomization
      print("> Residual Randomization")
      # Dantzig M
      S <- (t(x) %*% x) / n
      lambda_min <- min(.1, 2 * sqrt(log(p) / n))
      dantzig_M <- fastclime(S, lambda.min=lambda_min, nlambda=abs(nsolve))
      M = dantzig_M$icovlist[[length(dantzig_M$icovlist)]]
      out_rr = riSingle_hdi_dantzig(y, x, t(M), ind_0, beta_0, ind_1, beta_1, g_design)

      ci_rr_a = out_rr$ci_a
      ci_rr_n = out_rr$ci_n
      stopifnot(nrow(ci_rr_a)==length(beta_0))
      stopifnot(nrow(ci_rr_n)==length(beta_1))
      intv_rr_a = (ci_rr_a[,1] <= beta_0) & (ci_rr_a[,2] >= beta_0)
      intv_rr_n = (ci_rr_n[,1] <= beta_1) & (ci_rr_n[,2] >= beta_1)
      cov_rr_a = mean(intv_rr_a)
      cov_rr_n = mean(intv_rr_n)
      len_rr_a = mean(ci_rr_a[,2] - ci_rr_a[,1])
      len_rr_n = mean(ci_rr_n[,2] - ci_rr_n[,1])
      
      ############### Results ############### 
      if (nsolve < 0){
        Q = rbind(Q, c(cov_rr_a, len_rr_a,
                       cov_rr_n, len_rr_n))
      }
      else{
        Q = rbind(Q, c(cov_blpr_a, cov_hdi_a, cov_dlasso_a, cov_silm_a, cov_rr_a,
                       len_blpr_a, len_hdi_a, len_dlasso_a, len_silm_a, len_rr_a,
                       cov_blpr_n, cov_hdi_n, cov_dlasso_n, cov_silm_n, cov_rr_n,
                       len_blpr_n, len_hdi_n, len_dlasso_n, len_silm_n, len_rr_n))
      }
      print(Q)
      print(">> Coverage and Length")
      print(colMeans(Q))
    }
    
    ## Replications over. Save results.
    print("--===-    RESULTS ==--===")
    cm = as.data.frame(t(colMeans(Q, na.rm=T)))
    R1 = cbind(cm, data.frame(s=s0, x=X_design, b=beta_design, e=err_design, g=g_design, nsolve=nsolve))
    Results = rbind(Results, R1)
    write.csv(Results, file=sprintf("out/nsim_%d_s_%d_n_%d_p_%d_%s_%s_%s_%s_%d.csv",
                                    nsim, s0, n, p, X_design, beta_design, err_design, g_design, nsolve))
    print(Results)
    print("--===- --=-=== ==--===")
  }
}

# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$p_design, opt$nsim, opt$n, opt$p, opt$nsolve)
# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$p_design, opt$nsim, c(50, 100), c(100, 300), opt$nsolve)

# Test
main_sim(3, 'N2', 'D1', 'HMG', 'perm', 5, 50, 100, 300)



