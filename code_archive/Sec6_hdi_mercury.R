## Randomization method for high-dimensional inference.
##
require(hdi)
require(glmnet)
Rcpp::sourceCpp("cRAMSES_V.cpp")
source("RAMSES_V.R") # Randomization inference code.
require(mvtnorm)

args <- commandArgs(TRUE)  
JOB_ID = as.numeric(args[1])  # JOB_ID in 1-1200   (12 designs, 100 sim/design)
ERROR_TYPE = args[2]

# JOB_ID = JOB_ID + 1 # indexing from zero?

if(!is.na(JOB_ID)) { set.seed(JOB_ID) }

translate_jobid = function(jid) {
  nsim = 50
  stopifnot(jid %in% seq(1, nsim * 12))
  design_id = sim_id = NA
  if(jid %% nsim == 0) {
    design_id = jid/nsim; 
    sim_id = jid/design_id
  } else {
    sim_id = (jid %% nsim)
    design_id = 1+(jid - sim_id) / nsim
  }
  list(design_id=design_id, sim_id=sim_id)
}

translate_designid = function(did) {
  stopifnot(did %in% 1:12) ## 12 designs in total.
  s0 = beta_design = NA
  all_beta = c("U(0,2)", "U(0,4)", "U(-2,2)", "1", "2", "10")
  if(did > 6) {
    s0 = 15
    beta_design = all_beta[did-6]
  } else {
    s0 = 3
    beta_design = all_beta[did]
  }
  return(list(s0=s0, beta_design=beta_design))
}


get_Sigma = function(p) {
  file  = sprintf("sigma_p%d.rda", p)
  if(!file.exists(file)) {
    Sigma = matrix(apply(expand.grid(1:p, 1:p), 1, function(row) .9^abs(diff(row))), nrow=p)
    save(Sigma, file=file)
    return(Sigma)
  } else {
    load(file); return(Sigma)
  }
}

is.active = function(u) abs(u) > 1e-12

err_mixture = function(nsampl) {
  # a = sample(c(-1, 1), size=nsampl, T)* sample(log(1:nsampl))
  # a/sd(a)
  B  =rbinom(nsampl, size=1, prob=0.5)
  a = B * rnorm(nsampl, mean=-5, sd=0.5) + (1-B) * rnorm(nsampl, 5, 0.5)
  a/5
}


riTest_hdi = function(y, x) {
  n = length(y); p = ncol(x)
  
  #LASSO, RIDGE estimates
  beta_lasso = coef(cv.glmnet(x, y, intercept=F))[-1]
  out = solve_ridge(y, x)
  beta_ridge = out$bhat
  e = y - x %*% beta_lasso
  
  tn = function(u) out$P_inv %*% t(x) %*% u
  
  Tvals = matrix(0, nrow=p)
  for(R in 1:1000) {
    s = sample(c(-1, 1), n, T)
    enew = s * sample(e)  # support heteroskedastic errors.
    Tvals = cbind(Tvals, tn(enew))
  }
  
  lam0 = 0
  # Estimate the nuisance term.
  #
  # n0 = out$lam* c(t(lam)%*% out$P_inv %*% beta0)
  # est = tn(e) - bhat + lam0
  # est2 = out$lam* out$P_inv %*% beta_lasso
  # est = out$lam* out$P_inv %*% beta0   ## true value --> great results.
  est =  out$lam* out$P_inv %*% beta_lasso
  # print(sprintf("nuisance=%.2f -- est = %.2f", n0, est))
  
  Tobs =  beta_ridge - lam0 + est
  pvals = sapply(1:p, function(j) out_pval(list(tobs=Tobs[j], tvals=Tvals[j, ]), T))
  rej = (pvals <= .025/length(pvals))
  list(rej=rej, pvals=pvals)
}

solve_lasso= function(y, x) {
  stopifnot(length(y) == nrow(x))
  n = length(y)
  x1 = scale(x,TRUE,TRUE)/sqrt(n-1)
  # 
  out = lars(x1, y, intercept=F)
  E  = apply(coef(out), 1, function(row) which(is.active(row)))
  out1 = cv.glmnet(x1, y, intercept=F)
  c1 = which(is.active(coef(out1))) -1
  sel = which(unlist(lapply(E, function(l) setequal(l, c1))))
  
  if(length(sel) == 0) { return(list()) }
  
  # sel = selected coefficients based on CV.
  # out = cv.lars(x, y)
  bhat = coef(out)[sel,]; 
  lam1 = out$lambda[sel]
  # return the best-CV LASSO solution.
  list(lambda=lam1, betahat=bhat, x1=x1, S=which(is.active(bhat)))
}


solve_ridge = function(y, x) {
  stopifnot(length(y)==nrow(x))
  out  = cv.glmnet(x, y, alpha=0, intercept=F)
  lam = out$lambda.min
  p = ncol(x)
  I = diag(p)
  P = t(x) %*% x + lam * I; P_inv = solve(P)
  #
  bhat = P_inv %*% t(x) %*% y
  list(bhat=bhat, lam=lam, P=P, P_inv=P_inv)
}

old_rtest = function(y, x, num_R=100, H0=0) {
  
  n = length(y)
  # 1. get LASSO solution
  fit = solve_lasso(y, x)
  bhat = fit$betahat
  lam1 = fit$lambda
  x1 = fit$x1
  S = fit$S # observed support.

  # 2. prepare data structures to define test statistics.
  beta_S = bhat[S]
  # sign_S = sign(beta_S) ## fixed
  pred = x1 %*% bhat  # y^
  Xs = x1[, S]
  A_inv = solve(t(Xs) %*% Xs)
  C = A_inv %*% t(Xs) %*% x[, S]; C_inv = solve(C)
  D = C_inv %*% A_inv %*% t(Xs)
  
  ## Verify.
  ##
  # beta_S  # should be equal to....
  # as.numeric(A_inv %*% t(Xs) %*% y - lam1 * A_inv %*% sign_S)
  # as.numeric(C %*% beta0[S] + A_inv %*% t(Xs) %*% eps - lam1 * A_inv %*% sign_S)
  # # and the following should be equal.
  # C_inv %*%(beta_S + lam1 * A_inv %*% sign_S)
  # beta0[S] + C_inv %*% A_inv %*% t(Xs) %*% eps
  ##
  ## residuals
  er = y - pred
  print("Selected vars.")
  print(S)
  ## 3. Test all coefficients.
  pvals = c()
  for(j in 1:length(beta_S)) {
    # Test H0: beta_{S,j} = 0, beta_{-S} = 0.
    #
    # print(sprintf("Testing var. %d, true = %.2f, H0=%.3f", j, beta0[S][j], H0))
    proj = rep(0, length(beta_S)); proj[j] = 1
    proj = matrix(proj, nrow=1)  # = (0, 0, ...1, 0...0)
    
    tn = function(u) proj %*% (D %*% u)  # test statistic
    # observed test statistic
    testStat = function(u) {
      stopifnot(length(u)==length(S))
      as.numeric(proj %*% C_inv %*% (u + lam1 * A_inv %*% sign(u)) - H0)
    }
    ## observed test statistic
    Tobs = testStat(beta_S)
    # stopifnot(tn(eps) == tobs)  # this is true if null is true.
    print(Tobs)
    tvals = c()
    t0 = proc.time()[3]
    for(R in 1:num_R) {
      g = sample(c(-1, 1), size=n, replace=T)
      er_new = g * sample(er)
      y_new = pred + er_new
      out_new = solve_lasso(y_new, x)
      if(length(out_new) > 0) {
        # out = cv.lars(x, y)
        S_new = out_new$S
        
        if(setequal(S, S_new)) {
          tvals = c(tvals, tn(er_new))
        } else {
          # print(R)
          # print(S_new); print(S)
        }
        t1 = proc.time()[3]
        if(t1 - t0  > 5) {
          t0= t1
          print(sprintf("R=%d/%d--for var %d, util=%.2f%%...pval=%.4f", R, num_R, S[j], 100 * length(tvals)/R, 
                        out_pval(list(tobs=Tobs, tvals=tvals), T)))
          if(length(tvals) > 0) {
            hist(tvals, breaks=20)
            abline(v=Tobs, col="red", lwd=3)
            print(Tobs)
          }
        }
      }
    }
    
    print(sprintf("for var %d, util=%.2f%%...pval=%.4f", j, 100 * length(tvals)/num_R, 
                  out_pval(list(tobs=Tobs, tvals=tvals), T)))
    
    pvals = c(pvals, out_pval(list(tobs=Tobs, tvals=tvals), T))
  }
  
  ## 
  decision = (pvals <= (0.05/length(pvals)))
  return(S[which(decision)])
}


main_sim = function(s0, beta_design, job_id, n=100, p=500) {
  cols = c( "pwr_lasso", "pwr_rand", "fwer_lasso", "fwer_rand", "s0", "betas")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  Results = as.data.frame(Results)
  t0 = proc.time()[3]
  
  Sigma = get_Sigma(p)
  ## s0 = cardinality
  # x = rmvnorm(n, mean=rep(0, p), sigma = Sigma)
  # Z = matrix()
  # Sample design. Toeplitz
  x = rmvnorm(n, mean=rep(0, p), sigma = Sigma)
  Z = matrix()
  
  # change the support.
  S0 = sample(1:p, size=s0, replace=F)  # true support.
  Snot = setdiff(1:p, S0)

  beta = NA
  if(beta_design=="U(0,2)") {
    beta =  runif(p, min=0, max=2); 
  } else if(beta_design=="U(0,4)") {
    beta =  runif(p, min=0, max=4); 
  } else if(beta_design=="U(-2,2)") {
    beta = runif(p, min=-2, max=2)
  } else if(beta_design=="1") {
    beta = rep(1, p)
  } else if(beta_design=="2") {
    beta = rep(2, p)
  } else if(beta_design=="10") {
    beta = rep(10, p)
  }
  # knock out inactive set.
  beta[-S0] = 0
  
  ### Now is the main Iteration.
  M = matrix(0, nrow=0, ncol=4)  # interim results.
  colnames(M) = colnames(Results)[1:4]
  ###
  for(irep in 1:100) {
    ##
    ##
    eps = NA
    if(ERROR_TYPE=="normal") {
      eps = rnorm(n)
    } else if(ERROR_TYPE=="cauchy") {
      eps = rcauchy(n)
    } else if(ERROR_TYPE=="mixture") {
      eps = err_mixture(n)
    } else if(ERROR_TYPE=="t") {
      eps = rt(n, 3)
    }
    else {
      stop("Not supported error type.")
    }
    
    y = x %*% beta + eps
    
    if(nrow(Z) < 2) {
      lasso = lasso.proj(x, y, return.Z=T)
      Z = lasso$Z
      cache = list(x=x, z=Z)
      uid = as.integer(1e5*runif(1))
      save(cache, file=sprintf("out/cache/hdi_cache%d.rda", uid))
    }
    
    ## lasso solution
    lasso = lasso.proj(x, y, Z = Z)
    lasso_sel = which(lasso$pval.corr < 0.05) # LASSO selections.
    pwr_lasso = length(intersect(lasso_sel, S0)) / s0
    fwer_lasso = length(intersect(lasso_sel, Snot)) / (p-s0)
    
    ## Rand-test
    rand = riTest_hdi(y, x) ## vector of 1, 0, ...
    # print(rand$pvals[S0])
    pwr_rand = length(intersect(which(rand$rej==1), S0)) /s0
    fwer_rand = length(intersect(which(rand$rej==1), Snot))/(p - s0)
    
    # interim results.
    M = rbind(M, c(pwr_lasso, pwr_rand, fwer_lasso, fwer_rand))
    # Results = rbind(Results, c(pwr_lasso, pwr_rand, fwer_lasso, fwer_rand, s0, beta_design, ht))
    
    t1 = proc.time()[3]
    if(t1 - t0 > 10) {
      print(sprintf("iteration %d/%d, s0=%d, beta=%s, error-type = %s", irep, 100, s0, beta_design, ERROR_TYPE))
      print(M)
      print(colMeans(M))
      # Store results.
      t0 = t1
    }
  }
  ## Replications over. Save results.
  print("--===-    RESULTS ==--===")
  cm = as.data.frame(t(colMeans(M, na.rm=T)))
  R1 = cbind(cm, data.frame(s0=s0, beta_design=beta_design, error=ERROR_TYPE))
  Results = rbind(Results, R1)
  save(Results, file=sprintf("out/sim_hdi_reps%d_p%d_n%d__JOBID%d.rda", 100, p, n, job_id))
  print(Results)
  print("--===- --=-=== ==--===")
  
}

## Main simulation
if(!is.na(JOB_ID)) {
  ids = translate_jobid(JOB_ID)
  design = translate_designid(ids$design_id)
  ### main simulation.
  main_sim(design$s0, design$beta_design, JOB_ID)
}
