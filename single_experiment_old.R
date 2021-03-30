## RRI simulations
rm(list=ls())
pacman::p_load(fastclime, glmnet, lcmix, mvtnorm, optparse, parallel)

pacman::p_load(HDCI) # Liu et. al
pacman::p_load(hdi) # Buhlmann
# source("lasso_inference.R") # Montanari
# source("silm.R") # Zhang and Cheng
source("randomization_hdi.R") # RR

option_list = list(
  make_option(c("-s", "--sparsity"), type="integer"),
  make_option(c("-x", "--x_design"), type="character"),
  make_option(c("-b", "--b_design"), type="character"),
  make_option(c("-e", "--e_design"), type="character"),
  make_option(c("-g", "--g_design"), type="character"),
  make_option(c("-r", "--nsim"), type="integer"),
  make_option(c("-c", "--seed_index"), type="integer"),
  make_option(c("-n", "--n"), type="integer"),
  make_option(c("-p", "--p"), type="integer"),
  make_option(c("-d", "--n_draws"), type="integer"),
  make_option(c("-i", "--n_solve"), type="integer"),
  make_option(c("-m", "--scale"), type="integer")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opt = parse_args(opt_parser)


main_sim = function(s0, X_design, beta_design, err_design, g_design, nsim, seed_index, n_list, p_list, n_draws, n_solve, scale=TRUE){
  for (i in 1:length(n_list)){
    n = n_list[i]
    p = p_list[i]
    # p = 2 * n
    
    print(sprintf("--===-  s:           %s ==--===", s0))
    print(sprintf("--===-  X_design:    %s ==--===", X_design))
    print(sprintf("--===-  beta_design: %s ==--===", beta_design))
    print(sprintf("--===-  err_design:  %s ==--===", err_design))
    print(sprintf("--===-  g_design:    %s ==--===", g_design))
    print(sprintf("--===-  nsim:        %s ==--===", nsim))
    print(sprintf("--===-  n:           %s ==--===", n))
    print(sprintf("--===-  p:           %s ==--===", p))
    print(sprintf("--===-  n_draws:     %s ==--===", n_draws))
    print(sprintf("--===-  n_solve:     %s ==--===", n_solve))
    print(sprintf("--===-  scale:       %s ==--===", scale))
    
    n_cores = 1 
    
    if (n_solve < 0){
      # Name columns
      # cols = c( "a_cov_rr_10", "a_cov_rr_100", "a_cov_rr_250", "a_cov_rr_500", "a_cov_rr_750", "a_cov_rr_1000", "a_cov_rr_2000", "a_cov_rr_3000", "a_cov_rr_4000", "a_cov_rr_5000", "a_cov_rr_7500", "a_cov_rr_10000", 
      #           "a_len_rr_10", "a_len_rr_100", "a_len_rr_250", "a_len_rr_500", "a_len_rr_750", "a_len_rr_1000", "a_len_rr_2000", "a_len_rr_3000", "a_len_rr_4000", "a_len_rr_5000", "a_len_rr_7500", "a_len_rr_10000",
      #           "n_cov_rr_10", "n_cov_rr_100", "n_cov_rr_250", "n_cov_rr_500", "n_cov_rr_750", "n_cov_rr_1000", "n_cov_rr_2000", "n_cov_rr_3000", "n_cov_rr_4000", "n_cov_rr_5000", "n_cov_rr_7500", "n_cov_rr_10000",
      #           "n_len_rr_10", "n_len_rr_100", "n_len_rr_250", "n_len_rr_500", "n_len_rr_750", "n_len_rr_1000", "n_len_rr_2000", "n_len_rr_3000", "n_len_rr_4000", "n_len_rr_5000", "n_len_rr_7500", "n_len_rr_10000")
      cols = c( "a_cov_rr_10000", 
                "a_len_rr_10000",
                "n_cov_rr_10000",
                "n_len_rr_10000")
      # Cache results
      Q = matrix(0, nrow=0, ncol=length(cols)) 
      colnames(Q) = cols
    }
    else{
      cols = c("a_cov_blpr", "a_cov_dlasso", "a_cov_silm", "a_cov_rr",
               "a_len_blpr", "a_len_dlasso", "a_len_silm", "a_len_rr",
               "n_cov_blpr", "n_cov_dlasso", "n_cov_silm", "n_cov_rr", 
               "n_len_blpr", "n_len_dlasso", "n_len_silm", "n_len_rr")
      Q = matrix(0, nrow=0, ncol=length(cols)) 
      colnames(Q) = cols
    }
    Results = matrix(0, nrow=0, ncol=length(cols))
    colnames(Results) = cols
    Results = as.data.frame(Results)
    
    rho = 0.8
    Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
    
    for (i in (seed_index+1):(seed_index+nsim)){
      set.seed(i)
      
      print(sprintf("--===-  %d/%d iter ==--===",i, (seed_index+nsim)))
      ######################### START DGP ############################
      ## Generate X
      if(X_design=="N1"){
        x = mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=diag(p))
      } else if(X_design=="G1"){
        x = rmvgamma(n, shape=1, rate=1, corr=diag(p)) - 1
      } else if(X_design=="N2"){
        x = matrix((sample(c(-2, 2), size=n*p, replace=T) + rnorm(n*p)), nrow=n, byrow=TRUE)
      } else if(X_design=="TG"){
        x = mvtnorm::rmvnorm(n, mean=rep(0, p), sigma = Tptz)
      } else if(X_design=="TGM"){
        x = rmvgamma(n, shape=1, rate=1, corr=Tptz) - 1
      } else if(X_design=='L1'){
        x = rmvl(n, mu=rep(0, p), Sigma=diag(p))
      } else if(X_design=='TL'){
        x = rmvl(n, mu=rep(0, p), Sigma=Tptz)
      } else if(X_design=='T1'){
        x = rmvt(n, mu=rep(0, p), Sigma=diag(p))
      } else if (X_design=='TT'){
        x = rmvt(n, mu=rep(0, p), Sigma=Tptz)
      } else if(X_design=='WB'){
        x = rmvweisd(n, shape=rep(0.5, p), decay=1, corr=diag(p)) - gamma(2)
      # } else if(X_design=='TWB'){
      #   shape = rep(1, p)
      #   x = rmvweisd(n, shape=0.5, decay=1, corr=Tptz) - gamma(2)
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
      # ind_1 = ind[(s0+1):(s0+1)]
      beta_1 = beta[ind_1]
      
      ## Generate errors.
      if(err_design=="N1"){
        err = rnorm(n)
      } else if (err_design=="G1"){
        err = rgamma(n, shape=1, rate=1) - 1
      } else if(err_design=="N2"){
        err = sample(c(-2, 2), size=n, replace=T) + rnorm(n)
      } else if(err_design=="HG"){
        err = mvtnorm::rmvnorm(n=1, mean=rep(0, n), sigma=2 * diag(rowSums(x*x)/p))
      } else if(err_design=="HMG"){
        err = sample(c(-2, 2), size = n, replace = T) + 
          mvtnorm::rmvnorm(n=1, mean=rep(0, n), sigma=2 * diag(rowSums(x*x)/p))
      } else if(X_design=='L1'){
        err = rmvl(n, mu=0, Sigma=1)
      } else if(X_design=='T1'){
        err = rmvt(n, mu=0, Sigma=1)
      } else if(err_design=='WB'){
        err = rweisd(n, shape=0.5, 1) - gamma(2)
      }
      
      y = x %*% beta + as.vector(err)
      ######################### END DGP ############################
      
      ##################### START EXPERIMENT #######################
      # Residual Randomization
      print("> Residual Randomization")
      S <- t(x) %*% (x) / n
      # Solve for M
      lambda <- 0.25 * sqrt(log(p) / n)
      # Get path of Ms from fastclime
      clime_M <- fastclime(S, lambda.min=lambda, nlambda=abs(n_solve))
      
      out_rr = rr_min_clime(y, x, n_draws, clime_M, lambda, ind_0, beta_0, ind_1, beta_1, g_design, scale)
      
      # ci_rr_a = out_rr$ci_a
      # ci_rr_n = out_rr$ci_n
      # stopifnot(nrow(ci_rr_a)==length(beta_0))
      # stopifnot(nrow(ci_rr_n)==length(beta_1))
      # intv_rr_a = (ci_rr_a[,1] <= beta_0) & (ci_rr_a[,2] >= beta_0)
      # intv_rr_n = (ci_rr_n[,1] <= beta_1) & (ci_rr_n[,2] >= beta_1)
      # cov_rr_a = mean(intv_rr_a)
      # cov_rr_n = mean(intv_rr_n)
      # len_rr_a = mean(ci_rr_a[,2] - ci_rr_a[,1])
      # len_rr_n = mean(ci_rr_n[,2] - ci_rr_n[,1])
    
      cov_rr_a <- list()
      cov_rr_n <- list()
      len_rr_a <- list()
      len_rr_n <- list()
      
      for (i in names(out_rr$ci_a)){
        ci_rr_a = out_rr$ci_a[[i]]
        ci_rr_n = out_rr$ci_n[[i]]
        stopifnot(nrow(ci_rr_a)==length(beta_0))
        stopifnot(nrow(ci_rr_n)==length(beta_1))
        intv_rr_a = (ci_rr_a[,1] <= beta_0) & (ci_rr_a[,2] >= beta_0)
        intv_rr_n = (ci_rr_n[,1] <= beta_1) & (ci_rr_n[,2] >= beta_1)
        cov_rr_a[[i]] = mean(intv_rr_a)
        cov_rr_n[[i]] = mean(intv_rr_n)
        len_rr_a[[i]] = mean(ci_rr_a[,2] - ci_rr_a[,1])
        len_rr_n[[i]] = mean(ci_rr_n[,2] - ci_rr_n[,1])
      }
      
      rm(out_rr, clime_M)
      gc()
      
      # print(sprintf("1Norm(Beta): %f, 1Norm(dBeta): %f, Supp(Beta): %s, 1Norm(eps): %f",
      #               out_rr$norm1_beta, out_rr$norm1_dbeta, out_rr$norm0_beta, out_rr$norm1_eps))
      
      if (n_solve > 0) {
        ## 1. Residual Bootstrap Lasso + Partial Ridge
        print("> Res. Bootstrap + Ridge")
        out_blpr = bootLPR(x = x, y = y, type.boot = "residual", B = n_draws,
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

        rm(out_blpr)
        gc()

        # ## 2. Buhlmann, have experienced a weird crash, see slurm-1891414.out
        # print("> Bulhmann")
        # # Use wild when g_design is sign
        # if (g_design == "sign") {
        #   wild = TRUE
        # }
        # else{
        #   wild = FALSE
        # }
        # out_hdi = boot.lasso.proj(x, y, family = "gaussian", standardize = TRUE,
        #                           multiplecorr.method = "WY",
        #                           parallel = TRUE, ncores = n_cores,
        #                           betainit = "cv lasso", sigma = NULL, Z = NULL, verbose = FALSE,
        #                           return.Z = FALSE, robust= FALSE,
        #                           B = n_draws, boot.shortcut = FALSE,
        #                           return.bootdist = TRUE, wild = wild,
        #                           gaussian.stub = FALSE)
        # 
        # ci_hdi_a = confint(out_hdi, level = 0.95, parm=ind_0)
        # ci_hdi_n = confint(out_hdi, level = 0.95, parm=ind_1)
        # stopifnot(nrow(ci_hdi_a)==length(beta_0))
        # stopifnot(nrow(ci_hdi_n)==length(beta_1))
        # intv_hdi_a = (ci_hdi_a[,1] <= beta_0) & (ci_hdi_a[,2] >= beta_0)
        # intv_hdi_n = (ci_hdi_n[,1] <= beta_1) & (ci_hdi_n[,2] >= beta_1)
        # cov_hdi_a = mean(intv_hdi_a)
        # cov_hdi_n = mean(intv_hdi_n)
        # len_hdi_a = mean(ci_hdi_a[,2] - ci_hdi_a[,1])
        # len_hdi_n = mean(ci_hdi_n[,2] - ci_hdi_n[,1])
        # 
        # rm(out_hdi)
        # gc()

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

        rm(out_dlasso)
        gc()
        
        ## 4. SILM
        print("> SILM")
        ci_silm = SimE.CI(x, y, ind_0, ind_1, M = n_draws, alpha = 0.95)
        
        stopifnot(nrow(ci_silm$a_band.nst)==length(beta_0))
        stopifnot(nrow(ci_silm$n_band.nst)==length(beta_1))
        intv_silm_a = (ci_silm$a_band.nst[,1] <= beta_0) & (ci_silm$a_band.nst[,2] >= beta_0)
        intv_silm_n = (ci_silm$n_band.nst[,1] <= beta_1) & (ci_silm$n_band.nst[,2] >= beta_1)
        cov_silm_a = mean(intv_silm_a)
        cov_silm_n = mean(intv_silm_n)
        len_silm_a = mean(ci_silm$a_band.nst[,2] - ci_silm$a_band.nst[,1])
        len_silm_n = mean(ci_silm$n_band.nst[,2] - ci_silm$n_band.nst[,1])
        
        rm(ci_silm)
        gc()
      }
      
      ############### Results ############### 
      if (n_solve < 0){
        # Q = rbind(Q, c(cov_rr_a[['10']], cov_rr_a[['100']], cov_rr_a[['250']], cov_rr_a[['500']], cov_rr_a[['750']], cov_rr_a[['1000']], cov_rr_a[['2000']], cov_rr_a[['3000']],  cov_rr_a[['4000']],  cov_rr_a[['5000']],  cov_rr_a[['7500']],  cov_rr_a[['10000']], 
        #                len_rr_a[['10']], len_rr_a[['100']], len_rr_a[['250']], len_rr_a[['500']], len_rr_a[['750']], len_rr_a[['1000']], len_rr_a[['2000']], len_rr_a[['3000']],  len_rr_a[['4000']],  len_rr_a[['5000']],  len_rr_a[['7500']],  len_rr_a[['10000']], 
        #                cov_rr_n[['10']], cov_rr_n[['100']], cov_rr_n[['250']], cov_rr_n[['500']], cov_rr_n[['750']], cov_rr_n[['1000']], cov_rr_n[['2000']], cov_rr_n[['3000']],  cov_rr_n[['4000']],  cov_rr_n[['5000']],  cov_rr_n[['7500']],  cov_rr_n[['10000']], 
        #                len_rr_n[['10']], len_rr_n[['100']], len_rr_n[['250']], len_rr_n[['500']], len_rr_n[['750']], len_rr_n[['1000']], len_rr_n[['2000']], len_rr_n[['3000']],  len_rr_n[['4000']],  len_rr_n[['5000']],  len_rr_n[['7500']],  len_rr_n[['10000']]))
      }
         Q = rbind(Q, c(cov_rr_a[['10000']], 
                       len_rr_a[['10000']], 
                       cov_rr_n[['10000']], 
                       len_rr_n[['10000']]))
      }
      else{
        Q = rbind(Q, c(cov_blpr_a, cov_dlasso_a, cov_silm_a, cov_rr_a,
                       len_blpr_a, len_dlasso_a, len_silm_a, len_rr_a,
                       cov_blpr_n, cov_dlasso_n, cov_silm_n, cov_rr_n,
                       len_blpr_n, len_dlasso_n, len_silm_n, len_rr_n))
      }
      # print(Q)
      print(">> Coverage and Length")
      print(colMeans(Q))
    }
    
    ## Save and print results.
    print("--===-    RESULTS ==--===")
    # cm = as.data.frame(t(colMeans(Q, na.rm=T)))
    cm = as.data.frame(Q)
    R1 = cbind(cm, data.frame(s=s0, x=X_design, b=beta_design, e=err_design, g=g_design, n_draws=n_draws, n=n, p=p, n_solve=n_solve, scale=scale))
    Results = rbind(Results, R1)
    dir.create('out', showWarnings = FALSE)
    path  <- sprintf("out/rr_s_%d_x_%s_b_%s_e_%s_g_%s_r_%d_n_%d_p_%d_d_%d_i_%d_scale_%s", 
                     s0, X_design, beta_design, err_design, g_design, nsim, n, p, n_draws, n_solve, scale)
    dir.create(path, showWarnings = FALSE)
    write.csv(Results, file.path(path, sprintf("i_%d.csv", (seed_index+1))))
    print(Results)
    print("--===- --=-=== ==--===")
  }
}

# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, c(opt$n), c(opt$p), opt$n_draws, opt$n_solve, opt$scale)
main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index, c(50), c(100), opt$n_draws, opt$n_solve, opt$scale)
main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index, c(100), c(300), opt$n_draws, opt$n_solve, opt$scale)

# # Test
# ptm <- proc.time()
# main_sim(3, 'WB', 'D1', 'WB', 'perm', 5, 2, c(50), c(100), 1000, -500, 1)
# tt <- proc.time() - ptm
