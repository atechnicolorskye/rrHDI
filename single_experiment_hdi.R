## RRI simulations
rm(list=ls())
pacman::p_load(fastclime, glmnet, LaplacesDemon, lcmix, mvtnorm, optparse, parallel)

pacman::p_load(hdi) # Buhlmann

option_list = list(
  make_option(c("-s", "--sparsity"), type="integer"),
  make_option(c("-x", "--x_design"), type="character"),
  make_option(c("-b", "--b_design"), type="character"),
  make_option(c("-e", "--e_design"), type="character"),
  make_option(c("-g", "--g_design"), type="character"),
  make_option(c("-r", "--nsim"), type="integer"),
  make_option(c("-i", "--seed_index"), type="integer"),
  make_option(c("-n", "--n"), type="integer"),
  make_option(c("-p", "--p"), type="integer"),
  make_option(c("-d", "--n_draws"), type="integer")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opt = parse_args(opt_parser)

main_sim = function(s0, X_design, beta_design, err_design, g_design, nsim, seed_index, n_list, p_list, n_draws){
  for (i in 1:length(n_list)){
    n = n_list[i]
    p = p_list[i]
    # p = 2 * n
    
    print(sprintf("--===-  s:           %d ==--===", s0))
    print(sprintf("--===-  X_design:    %s ==--===", X_design))
    print(sprintf("--===-  beta_design: %s ==--===", beta_design))
    print(sprintf("--===-  err_design:  %s ==--===", err_design))
    print(sprintf("--===-  g_design:    %s ==--===", g_design))
    print(sprintf("--===-  nsim:        %s ==--===", nsim))
    print(sprintf("--===-  seed_index:  %s ==--===", seed_index))
    print(sprintf("--===-  n:           %s ==--===", n))
    print(sprintf("--===-  p:           %s ==--===", p))
    print(sprintf("--===-  n_draws:     %s ==--===", n_draws))
    
    cols = c("a_cov_hdi", "a_len_hdi", "n_cov_hdi", "n_len_hdi")
    Q = matrix(0, nrow=0, ncol=length(cols)) 
    colnames(Q) = cols
    Results = matrix(0, nrow=0, ncol=length(cols))
    colnames(Results) = cols
    Results = as.data.frame(Results)
    
    rho = 0.8
    Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
    
    n_cores = 3
    
    for (i in (seed_index+1):(seed_index+nsim)){
      set.seed(i)
      
      print(sprintf("--===-  %d/%d iter ==--===", i, (seed_index+nsim)))
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
        shape = rep(1, p)
        x = rmvweisd(n, shape=shape, decay=shape) - gamma(2)
      } else if(X_design=='TWB'){
        shape = rep(1, p)
        x = rmvweisd(n, shape=shape, decay=shape, corr=Tptz) - gamma(2)
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
        err = rweisd(n, 1, 1) - gamma(2)
      }
      
      y = x %*% beta + as.vector(err)
      ######################### END DGP ############################
      
      ##################### START EXPERIMENT #######################
      ## 2. HDI, have experienced a weird crash, see slurm-1891414.out
      print("> HDI")
      # Use wild when g_design is sign
      if (g_design == "sign") {
        wild = TRUE
      }
      else{
        wild = FALSE
      }
      out_hdi = boot.lasso.proj(x, y, family = "gaussian", standardize = TRUE,
                                multiplecorr.method = "WY",
                                parallel = TRUE, ncores = n_cores,
                                betainit = "cv lasso", sigma = NULL, Z = NULL, verbose = FALSE,
                                return.Z = FALSE, robust= FALSE,
                                B = n_draws, boot.shortcut = FALSE,
                                return.bootdist = TRUE, wild = wild,
                                gaussian.stub = FALSE)
      
      ci_hdi_a = confint(out_hdi, level = 0.95, parm=ind_0)
      ci_hdi_n = confint(out_hdi, level = 0.95, parm=ind_1)
      stopifnot(nrow(ci_hdi_a)==length(beta_0))
      stopifnot(nrow(ci_hdi_n)==length(beta_1))
      intv_hdi_a = (ci_hdi_a[,1] <= beta_0) & (ci_hdi_a[,2] >= beta_0)
      intv_hdi_n = (ci_hdi_n[,1] <= beta_1) & (ci_hdi_n[,2] >= beta_1)
      cov_hdi_a = mean(intv_hdi_a)
      cov_hdi_n = mean(intv_hdi_n)
      len_hdi_a = mean(ci_hdi_a[,2] - ci_hdi_a[,1])
      len_hdi_n = mean(ci_hdi_n[,2] - ci_hdi_n[,1])
      
      rm(out_hdi)
      gc()
      
      ############### Results ############### 
      Q = rbind(Q, c(cov_hdi_a, len_hdi_a, cov_hdi_n, len_hdi_n))
      print(">> Coverage and Length")
      print(colMeans(Q))
    }
    
    ## Save and print results.
    print("--===-    RESULTS ==--===")
    cm = as.data.frame(Q)
    R1 = cbind(cm, data.frame(s=s0, x=X_design, b=beta_design, e=err_design, g=g_design, n_draws=n_draws, seed_index=(seed_index+1), n=n, p=p))
    Results = rbind(Results, R1)
    dir.create('out', showWarnings = FALSE)
    path  <- sprintf("out/hdi_s_%d_x_%s_b_%s_e_%s_g_%s_r_%d_n_%d_p_%d_d_%d", 
                     s0, X_design, beta_design, err_design, g_design, nsim, n, p, n_draws)
    dir.create(path, showWarnings = FALSE)
    write.csv(Results, file.path(path, sprintf("i_%d.csv", (seed_index+1))))
    print(Results)
    print("--===- --=-=== ==--===")
  }
}

# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index,, c(opt$n), c(opt$p), opt$n_draws)
# main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index, c(50), c(100), opt$n_draws)
main_sim(opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index, c(100), c(300), opt$n_draws)

# # Test
# ptm <- proc.time()
# main_sim(10, 'TG', 'D1', 'G1', 'perm', 2, 250, c(50), c(100), 1000)
# tt <- proc.time() - ptm
