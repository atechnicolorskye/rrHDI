## RRI simulations
rm(list=ls())
pacman::p_load(fastclime, glmnet, lcmix, mvtnorm, optparse, parallel)

pacman::p_load(HDCI) # Liu et. al
pacman::p_load(hdi) # Buhlmann
source("lasso_inference.R") # Montanari
source("silm.R") # Zhang and Cheng
source("randomization_hdi.R") # RR

option_list = list(
  make_option(c("-m", "--procedure"), type="character"),
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
  make_option(c("-v", "--standard"), type="integer")
); 

opt_parser = OptionParser(option_list=option_list, add_help_option=FALSE)
opt = parse_args(opt_parser)

main_sim = function(procedure, s0, X_design, beta_design, err_design, g_design, nsim, seed_index, n_list, p_list, n_draws, standard){
    for (i in 1:length(n_list)){
        n = n_list[i]
        p = p_list[i]
        
        print(sprintf("--===-  procedure:   %s ==--===", procedure))
        print(sprintf("--===-  s:           %s ==--===", s0))
        print(sprintf("--===-  X_design:    %s ==--===", X_design))
        print(sprintf("--===-  beta_design: %s ==--===", beta_design))
        print(sprintf("--===-  err_design:  %s ==--===", err_design))
        print(sprintf("--===-  g_design:    %s ==--===", g_design))
        print(sprintf("--===-  nsim:        %s ==--===", nsim))
        print(sprintf("--===-  seed_index:  %s ==--===", seed_index))
        print(sprintf("--===-  n:           %s ==--===", n))
        print(sprintf("--===-  p:           %s ==--===", p))
        print(sprintf("--===-  n_draws:     %s ==--===", n_draws))
        print(sprintf("--===-  standard:    %s ==--===", standard))
        
        cols = c(paste("a_cov_i_", procedure, sep=""),
                 paste("a_cov_a_", procedure, sep=""),
                 paste("a_cov_s_", procedure, sep=""),
                 paste("a_len_i_", procedure, sep=""),
                 paste("a_len_a_", procedure, sep=""),
                 paste("a_len_s_", procedure, sep=""),
                 paste("n_cov_i_", procedure, sep=""),
                 paste("n_cov_a_", procedure, sep=""),
                 paste("n_cov_s_", procedure, sep=""),
                 paste("n_len_i_", procedure, sep=""),
                 paste("n_len_a_", procedure, sep=""),
                 paste("n_len_s_", procedure, sep="")
                 )
        Q = matrix(0, nrow=0, ncol=length(cols)) 
        colnames(Q) = cols
        Results = matrix(0, nrow=0, ncol=length(cols))
        colnames(Results) = cols
        Results = as.data.frame(Results)
        
        rho = 0.8
        Tptz = rho^abs(matrix(1:p, nrow=p, ncol=p, byrow=F) - matrix(1:p, nrow=p, ncol=p, byrow=T))
        
        if (procedure == 'hdi'){
          n_cores = 3
        } else {
          n_cores = 1
        }
        
        for (i in (seed_index+1):(seed_index+nsim)){
            set.seed(i)
            
            print(sprintf("--===-  %d/%d iter ==--===",i, (seed_index+nsim)))
            ######################### START DGP ############################
            ## Generate X
            if (X_design=="N1"){
                x = mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=diag(p))
            } else if (X_design=="G1"){
                x = rmvgamma(n, shape=1, rate=1, corr=diag(p)) - 1
            } else if (X_design=="N2"){
                x = matrix((sample(c(-2, 2), size=n*p, replace=T) + rnorm(n*p)), nrow=n, byrow=TRUE)
            } else if (X_design=="TG"){
                x = mvtnorm::rmvnorm(n, mean=rep(0, p), sigma = Tptz)
            } else if (X_design=="TGM"){
                x = rmvgamma(n, shape=1, rate=1, corr=Tptz) - 1
            } else if (X_design=='L1'){
                x = rmvl(n, mu=rep(0, p), Sigma=diag(p))
            } else if (X_design=='TL'){
                x = rmvl(n, mu=rep(0, p), Sigma=Tptz)
            } else if (X_design=='T1'){
                x = rmvt(n, mu=rep(0, p), Sigma=diag(p))
            } else if (X_design=='TT'){
                x = rmvt(n, mu=rep(0, p), Sigma=Tptz)
            } else if (X_design=='WB'){
              x = rmvweisd(n, shape=rep(0.5, p), decay=1, corr=diag(p)) - gamma(2)
              # } else if(X_design=='TWB'){
              #   shape = rep(1, p)
              #   x = rmvweisd(n, shape=0.5, decay=1, corr=Tptz) - gamma(2)
            }
            
            ## Generate beta
            beta = rep(0, p)
            if (beta_design=="D1"){
                beta[1:s0] = sample(c(-1, 1), size=s0, replace=T)
            } else if (beta_design=="D5") {
                beta[1:s0] = sample(c(-5, 5), size=s0, replace=T)
            }
            # Shuffle
            beta[s0:(s0+1)] = beta[(s0+1):s0]
            ind = rev(order(abs(beta)))
            # Pick active variables, first isolated, second adjacent and third sandwiched
            ind_0 =  ind[1:3]
            beta_0 = beta[ind_0]
            # Pick 15 inactive variables, first adjacent and last two beside isolated
            ind_1 = ind[(p-2):p]
            beta_1 = beta[ind_1]
            
            ## Generate errors.
            if (err_design=="N1"){
                err = rnorm(n)
            } else if (err_design=="G1"){
                err = rgamma(n, shape=1, rate=1) - 1
            } else if (err_design=="N2"){
                err = sample(c(-2, 2), size=n, replace=T) + rnorm(n)
            } else if (err_design=="HG"){
                err = mvtnorm::rmvnorm(n=1, mean=rep(0, n), sigma=2 * diag(rowSums(x*x)/p))
            } else if (err_design=="HMG"){
                err = sample(c(-2, 2), size = n, replace = T) + 
                    mvtnorm::rmvnorm(n=1, mean=rep(0, n), sigma=2 * diag(rowSums(x*x)/p))
            } else if (X_design=='L1'){
                err = rmvl(n, mu=0, Sigma=1)
            } else if (X_design=='T1'){
                err = rmvt(n, mu=0, Sigma=1)
            } else if (err_design=='WB'){
                err = rweisd(n, 0.5, 1) - gamma(2)
            }
            
            y = x %*% beta + as.vector(err)
            
            if (standard == 1){
              sd <- attr(scale(x), "scaled:scale")
              x <- scale(x, center = FALSE)
            } else {
              sd <- rep(1, p)
            }
            
            ######################### END DGP ############################
            
            ##################### START EXPERIMENT #######################
            if (procedure=='rr'){
              print("> Residual Randomization")
              # S <- t(x) %*% (x) / n
              # # Solve for M
              # lambda <- 0.1 * sqrt(log(p) / n)
              # # Get path of Ms from fastclime
              # clime_M <- fastclime(S, lambda.min=lambda, nlambda=500)
              clime_M <- array(0, c(p, p))
              
              out_rr = rr_min_clime(y, x, n_draws, clime_M, lambda, ind_0, beta_0, ind_1, beta_1, g_design, 1)
              
              # cov_a <- list()
              # cov_n <- list()
              # len_a <- list()
              # len_n <- list()
              
              delta = '10000'
              # for (delta in c('1000', '2500', '5000', '7500', '10000')){
              ci_a = out_rr$ci_a[[delta]] / sd[ind_0] 
              ci_n = out_rr$ci_n[[delta]] / sd[ind_1]
              stopifnot(nrow(ci_a)==length(beta_0))
              stopifnot(nrow(ci_n)==length(beta_1))
              intv_a = (ci_a[,1] <= beta_0) & (ci_a[,2] >= beta_0)
              intv_n = (ci_n[,1] <= beta_1) & (ci_n[,2] >= beta_1)
              cov_a = intv_a
              cov_n = intv_n
              len_a = ci_a[,2] - ci_a[,1]
              len_n = ci_n[,2] - ci_n[,1]
                # }
              # rm(out_rr, clime_M)
              rm(out_rr)
              gc()
            
            } else if (procedure=='blpr'){
              print("> Res. Bootstrap + Ridge")
              out_blpr = bootLPR(x = x, y = y, type.boot = "residual", B = n_draws,
                                 parallel = TRUE, parallel.boot = TRUE, ncores.boot = n_cores)
  
              ci_a = t(out_blpr$interval.LPR)[ind_0,] / sd[ind_0]
              ci_n = t(out_blpr$interval.LPR)[ind_1,] / sd[ind_1]
              stopifnot(nrow(ci_a)==length(beta_0))
              stopifnot(nrow(ci_n)==length(beta_1))
              intv_a = (ci_a[,1] <= beta_0) & (ci_a[,2] >= beta_0)
              intv_n = (ci_n[,1] <= beta_1) & (ci_n[,2] >= beta_1)
              cov_a = intv_a
              cov_n = intv_n
              len_a = ci_a[,2] - ci_a[,1]
              len_n = ci_n[,2] - ci_n[,1]
  
              rm(out_blpr)
              gc()
              
            } else if (procedure=='hdi'){
              print("> Bulhmann")
              # Use wild when g_design is sign
              wild = FALSE
              if (g_design == "sign") {
                wild = TRUE
              }
              out_hdi = boot.lasso.proj(x, y, family = "gaussian", standardize = TRUE,
                                        multiplecorr.method = "WY",
                                        parallel = TRUE, ncores = n_cores,
                                        betainit = "cv lasso", sigma = NULL, Z = NULL, verbose = FALSE,
                                        return.Z = FALSE, robust= FALSE,
                                        B = n_draws, boot.shortcut = FALSE,
                                        return.bootdist = TRUE, wild = wild,
                                        gaussian.stub = FALSE)

              ci_a = confint(out_hdi, level = 0.95, parm=ind_0) / sd[ind_0]
              ci_n = confint(out_hdi, level = 0.95, parm=ind_1) / sd[ind_1]
              stopifnot(nrow(ci_a)==length(beta_0))
              stopifnot(nrow(ci_n)==length(beta_1))
              intv_a = (ci_a[,1] <= beta_0) & (ci_a[,2] >= beta_0)
              intv_n = (ci_n[,1] <= beta_1) & (ci_n[,2] >= beta_1)
              cov_a = intv_a
              cov_n = intv_n
              len_a = ci_a[,2] - ci_a[,1]
              len_n = ci_n[,2] - ci_n[,1]

              rm(out_hdi)
              gc()
              
            } else if (procedure=='dlasso'){
              print("> Debiased LASSO")
              out_dlasso = SSLasso(x, y, intercept = FALSE)
  
              ci_a = cbind(out_dlasso$low.lim[ind_0] / sd[ind_0], out_dlasso$up.lim[ind_0] / sd[ind_0])
              ci_n = cbind(out_dlasso$low.lim[ind_1] / sd[ind_1], out_dlasso$up.lim[ind_1] / sd[ind_1])
              stopifnot(nrow(ci_a)==length(beta_0))
              stopifnot(nrow(ci_n)==length(beta_1))
              intv_a = (ci_a[,1] <= beta_0) & (ci_a[,2] >= beta_0)
              intv_n = (ci_n[,1] <= beta_1) & (ci_n[,2] >= beta_1)
              cov_a = intv_a
              cov_n = intv_n
              len_a = ci_a[,2] - ci_a[,1]
              len_n = ci_n[,2] - ci_n[,1]
  
              rm(out_dlasso)
              gc()
              
            } else if (procedure=='silm'){
              print("> SILM")
              ci_silm = SimE.CI(x, y, ind_0, ind_1, M = n_draws, alpha = 0.95)
              
              ci_a = cbind(ci_silm$a_band.st[,1] / sd[ind_0], ci_silm$a_band.st[,2] / sd[ind_0])
              ci_n = cbind(ci_silm$n_band.st[,1] / sd[ind_1], ci_silm$n_band.st[,2] / sd[ind_1])
              stopifnot(nrow(ci_a)==length(beta_0))
              stopifnot(nrow(ci_n)==length(beta_1))
              intv_a = (ci_a[,1] <= beta_0) & (ci_a[,2] >= beta_0)
              intv_n = (ci_n[,1] <= beta_1) & (ci_n[,2] >= beta_1)
              cov_a = intv_a
              cov_n = intv_n
              len_a = ci_a[,2] - ci_a[,1]
              len_n = ci_n[,2] - ci_n[,1]
              
              # stopifnot(nrow(ci_silm$a_band.st)==length(beta_0))
              # stopifnot(nrow(ci_silm$n_band.st)==length(beta_1))
              # intv_silm_a = (ci_silm$a_band.st[,1] / sd[ind_0] <= beta_0) & (ci_silm$a_band.st[,2] / sd[ind_0] >= beta_0)
              # intv_silm_n = (ci_silm$n_band.st[,1] / sd[ind_1] <= beta_1) & (ci_silm$n_band.st[,2] / sd[ind_1] >= beta_1)
              # cov_silm_a = mean(intv_silm_a)
              # cov_silm_n = mean(intv_silm_n)
              # len_silm_a = mean(ci_silm$a_band.st[,2] / sd[ind_0] - ci_silm$a_band.st[,1] / sd[ind_0])
              # len_silm_n = mean(ci_silm$n_band.st[,2] / sd[ind_1] - ci_silm$n_band.st[,1] / sd[ind_1])
              
              rm(ci_silm)
              gc()
            }
            
            ############### Results ############### 
            Q = rbind(Q, c(cov_a,
                           len_a,
                           cov_n,
                           len_n)
                      )
            # print(Q)
            print(">> Coverage and Length")
            print(colMeans(Q))
        }
        
        ## Save and print results.
        print("--===-    RESULTS ==--===")
        # cm = as.data.frame(t(colMeans(Q, na.rm=T)))
        cm = as.data.frame(Q)
        R1 = cbind(cm, data.frame(s=s0, x=X_design, b=beta_design, e=err_design, g=g_design, n_draws=n_draws, n=n, p=p, standard=standard))
        Results = rbind(Results, R1)
        dir.create('out', showWarnings = FALSE)
        path  <- sprintf("out/%s_s_%d_x_%s_b_%s_e_%s_g_%s_r_%d_n_%d_p_%d_d_%d_v_%d", 
                         procedure, s0, X_design, beta_design, err_design, g_design, nsim, n, p, n_draws, standard)
        dir.create(path, showWarnings = FALSE)
        write.csv(Results, file.path(path, sprintf("i_%d.csv", (seed_index+1))))
        print(Results)
        print("--===- --=-=== ==--===")
    }
}

# main_sim(opt$procedure, opt$sparsity, opt$x_design, opt$b_design, opt$e_design, opt$g_design, opt$nsim, opt$seed_index, c(opt$n), c(opt$p), opt$n_draws, opt$standard)

# # Test
# ptm <- proc.time()
main_sim('rr', 4, 'WB', 'D1', 'N1', 'perm', 100, 0, c(50), c(100), 1000, 1)
# tt <- proc.time() - ptm
