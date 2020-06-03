## Reproduce results for Bernoulli 2013: Ridge Projection for hdi
rm(list=ls())
require(hdi)
library(glmnet)
library(mvtnorm)

main_sim = function(nsim=500, n=100){
  cols = c( "pwr_single", "fwer_single", "pwr_multi", "fwer_multi", "p", "Sigma", "s0", "betas")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  Results = as.data.frame(Results)
  t0 = proc.time()[3]
  
  all_Sigma_designs = c("id", "equi.corr")
  all_p_designs = c(500, 2500)
  all_s0_designs = c(15)
  # all_beta_designs = c("0.25", "0.5", "1")
  all_beta_designs = c("1")

    # p: number of covariates
    for(p in all_p_designs){
      # Sigma: covariance matrix ==> M1, M2
      for (Sigma in all_Sigma_designs){
        Sig = NA
        if(Sigma == "id"){
          Sig = diag(p)
        } else if(Sigma == "equi.corr"){
          Sig = .8 * matrix(1, p, p)
          diag(Sig) = 1
        }
        
        # s0: number of true support
        for (s0 in all_s0_designs) {
          S0 = sample(1:p, size=s0, replace=F)  # true support.
          Snot = setdiff(1:p, S0)
          # beta: coefficients
          for (beta_design in all_beta_designs) {
            beta = NA
            if(beta_design=="0.25"){
              beta = rep(0.25, p)
            } else if(beta_design=="0.5"){
              beta = rep(0.5, p)
            } else if(beta_design=="1"){
              beta = rep(1, p)
            }
            # knock out inactive set.
            beta[-S0] = 0
            
            ### Now is the main Iteration.
            M = matrix(0, nrow=0, ncol=4)  # interim results.
            colnames(M) = colnames(Results)[1:4]
            
            for (i in 1:nsim) {
              x = rmvnorm(n, mean=rep(0, p), sigma = Sig)
              eps = rnorm(n) # errors
              y = x %*% beta + eps
              
              ## ridge lasso solution
              # ridge = ridge.proj(x, y, betainit = "cv lasso", suppress.grouptesting = FALSE)
              ridge = ridge.proj(x, y, betainit = "scaled lasso", suppress.grouptesting = TRUE)
              
              ## Single testing results
              ridge_sel_single = which(ridge$pval < 0.05) # Single testing selections.
              pwr_ridge_single = length(intersect(ridge_sel_single, S0)) / s0
              fwer_ridge_single = length(intersect(ridge_sel_single, Snot)) / (p-s0)
              
              ## Multiple testing results
              ridge_sel_multi = which(ridge$pval.corr < 0.05) # Multiple testing selections.
              pwr_ridge_multi = length(intersect(ridge_sel_multi, S0)) / s0
              fwer_ridge_multi = length(intersect(ridge_sel_multi, Snot)) / (p-s0)
              
              # interim results.
              M = rbind(M, c(pwr_ridge_single, fwer_ridge_single, pwr_ridge_multi, fwer_ridge_multi))
              
              t1 = proc.time()[3]
              print(sprintf("Time %f", (t1-t0)))
              if(t1 - t0 > 10) {
                print(sprintf("iteration %d, s0=%d, beta=%s", i, s0, beta_design))
                print(M)
                print(colMeans(M))
                # Store results.
                t0 = t1
              }
            }
            ## Replications over. Save results.
            print("--===-    RESULTS ==--===")
            cm = as.data.frame(t(colMeans(M, na.rm=T)))
            R1 = cbind(cm, data.frame(p=p, Sig=Sigma, s0=s0, beta_design=beta_design))
            Results = rbind(Results, R1)
            save(Results, file=sprintf("out/sim_hdi_sims%d_n%d.rda", nsim, n))
            print(Results)
            print("--===- --=-=== ==--===")
        }
      }
    }
  }
}

main_sim()



