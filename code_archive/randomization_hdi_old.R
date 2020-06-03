library(glmnet)
library(hdi)
library(mvtnorm)
library(RPtests)

Rcpp::sourceCpp("cRAMSES_V.cpp")
source("RAMSES_V.R") # Randomization Inferenc

args <- commandArgs(TRUE)  # (H0, num_G, nreps=1000)

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

riCI_hdi = function(y, x, j, vals=NA) {
  if(is.na(vals)) {
    vals = seq(-1, 1, length.out=20)
  }
  for(v0 in vals) {
    stop("Implemenet")
  }
}

riSingle_hdi_new = function(y, x, M, S_test, lam0_vals, vals=NA) {
  # M = Montanari's matrix
  # Test H0_j: beta_j = lam0_j  simultaneously
  # S_test = (id1, id2...) ids of variables to test (to  speed up things.)
  #
  n = nrow(x); p=ncol(x)
  stopifnot(length(y)==nrow(x))

  s = length(S_test)
  Y = matrix(y, nrow=n, ncol=s, byrow=F) # make this nxp matrix
  C = matrix(0, nrow=s, ncol=p)  # s x p   binary indicator of  H0_j
  for(i in 1:s) { C[i, S_test[i]] = 1 }

  t0 = proc.time()[3]
  # test H0_j: beta_j = v0, for all j in S_test.
  # 1. Create restricted LASSO solutions.
  Beta_l = matrix(0, nrow=p, ncol=s)
  cnt = 0
  for(jj in 1:s) {
    j = S_test[jj] # jj = 1, 2...   j = index  from 1-p
    v0 = lam0_vals[jj]

    beta_lasso = coef(cv.glmnet(x[,-j], y - x[,j]*v0, intercept=F))[-1]
    # solve_lasso(y - x[,j]*v0, x[,-j])  # H0_j: beta_j = v0
    pre = c(); post = c();
    if(j > 1) { pre = beta_lasso[seq(1, j-1)] }
    if(j < p) { post = beta_lasso[seq(j, p-1)] }
    bhat_j = c(pre, v0, post)
    Beta_l[, jj] = bhat_j
  }

  # 2. restricted errors.
  Eps = Y - x %*% Beta_l  # n x s, Eps[, j] = errors under H0_j.
  An = (1/sqrt(n)) * C %*%  M %*% t(x)   # s x n.
  t_n = function(u) {
    diag(An %*% u)
  }
  Tobs = t_n(Eps)  # vector with s-elements.

  Tvals = matrix(0, nrow=0, ncol=s)   # R x s
  for(R in 1:999) {
    G = diag(sample(c(-1, 1), size=n, replace=T))
    Tvals = rbind(Tvals, t_n(G %*% Eps))
  }

  # 3. PVals.
  Pvals = sapply(1:s, function(j) {
    out_pval(list(tobs=Tobs[j], tvals=Tvals[,j]), ret_pval = T)
  })

  return(Pvals)
}


riCI_hdi_new = function(y, x, M, S_test, vals=NA) {
  # M = Montanari's matrix
  # Test H0_j: beta_j = a0_j  simultaneously
  # S_test = (id1, id2...) ids of variables to test (to  speed up things.)
  #
  n = nrow(x); p=ncol(x)
  stopifnot(length(y)==nrow(x))

  s = length(S_test)
  Y = matrix(y, nrow=n, ncol=s, byrow=F) # make this nxp matrix
  C = matrix(0, nrow=s, ncol=p)  # s x p   binary indicator of  H0_j
  for(i in 1:s) { C[i, S_test[i]] = 1 }

  test_vals = seq(-0.5, 1.2, length.out=20)
  Pvals = matrix(0, nrow=0, ncol=s)

  t0 = proc.time()[3]
  for(v0 in test_vals) {
    # test H0_j: beta_j = v0, for all j in S_test.
    # 1. Create restricted LASSO solutions.
    Beta_l = matrix(0, nrow=p, ncol=s)
    cnt = 0
    for(j in S_test) {
      cnt = cnt + 1
      beta_lasso = coef(cv.glmnet(x[,-j], y - x[,j]*v0, intercept=F))[-1]
      # solve_lasso(y - x[,j]*v0, x[,-j])  # H0_j: beta_j = v0
      pre = c(); post = c();
      if(j > 1) { pre = beta_lasso[seq(1, j-1)] }
      if(j < p) { post = beta_lasso[seq(j, p-1)] }
      bhat_j = c(pre, v0, post)
      Beta_l[,cnt] = bhat_j
      # print(cnt)
    }


    # beta_lasso = coef(cv.glmnet(x, y, intercept=F))[-1]
    # Beta_l = matrix(beta_lasso, nrow=p, ncol=s, byrow=F) # same LASSO solution

    # 2. restricted errors.
    Eps = Y - x %*% Beta_l  # n x s, Eps[, j] = errors under H0_j.
    An = (1/sqrt(n)) * C %*%  M %*% t(x)   # s x n.
    t_n = function(u) {
      diag(An %*% u)
    }
    Tobs = t_n(Eps)  # vector with s-elements.

    Tvals = matrix(0, nrow=0, ncol=s)   # R x s
    for(R in 1:1000) {
      G = diag(sample(c(-1, 1), size=n, replace=T))
      Tvals = rbind(Tvals, t_n(G %*% Eps))
    }

    # 3. PVals.
    pls = sapply(1:s, function(j) {
      out_pval(list(tobs=Tobs[j], tvals=Tvals[,j]), ret_pval = T)
    })
    # print(sprintf("Testing %.2f", v0))
    Pvals = rbind(Pvals, pls)
    t1 = proc.time()[3]
    if(t1 - t0 > 5) {
      print(sprintf("v = %.2f max=%.2f", v0, max(test_vals)))
      t0 = t1
    }
  }

  rownames(Pvals) = round(test_vals, 2)
  return(list(val=test_vals, pvals=Pvals))
}


riCI_hdi_classic  = function(y, x, S_test, vals=NA) {
  # M = Montanari's matrix
  # Test H0_j: beta_j = a0_j  simultaneously
  # S_test = (id1, id2...) ids of variables to test (to  speed up things.)
  #
  n = nrow(x); p=ncol(x)
  stopifnot(length(y)==nrow(x))

  s = length(S_test)
  Y = matrix(y, nrow=n, ncol=s, byrow=F) # make this nxp matrix
  C = matrix(0, nrow=s, ncol=p)  # s x p   binary indicator of  H0_j
  for(i in 1:s) { C[i, S_test[i]] = 1 }

  test_vals = seq(-0.5, 1.2, length.out=20)
  Pvals = matrix(0, nrow=0, ncol=s)

  t0 = proc.time()[3]
  for(v0 in test_vals) {
    # test H0_j: beta_j = v0, for all j in S_test.
    # 1. Create restricted LASSO solutions.
    Beta_l = matrix(0, nrow=p, ncol=s)
    cnt = 0
    for(j in S_test) {
      cnt = cnt + 1
      out_j = solve_lasso(y - x[,j]*v0, x[,-j])  # H0_j: beta_j = v0
      pre = c(); post = c();
      if(j > 1) { pre = out_j$betahat[seq(1, j-1)] }
      if(j < p) { post = out_j$betahat[seq(j, p-1)] }
      bhat_j = c(pre, v0, post)
      Beta_l[,cnt] = bhat_j
      # print(cnt)
    }
    # out = solve_lasso(y, x)  # H0_j: beta_j = v0
    # Beta_l = matrix(out$betahat, nrow=p, ncol=s, byrow=F) # same LASSO solution

    # 2. restricted errors.
    Eps = Y - x %*% Beta_l  # n x s, Eps[, j] = errors under H0_j.
    An = (1/sqrt(n)) * C %*%  M %*% t(x)   # s x n.
    t_n = function(u) {
      diag(An %*% u)
    }
    Tobs = t_n(Eps)  # vector with s-elements.

    Tvals = matrix(0, nrow=0, ncol=s)   # R x s
    for(R in 1:500) {
      G = diag(sample(c(-1, 1), size=n, replace=T))
      Tvals = rbind(Tvals, t_n(G %*% Eps))
    }

    # 3. PVals.
    pls = sapply(1:s, function(j) {
      out_pval(list(tobs=Tobs[j], tvals=Tvals[,j]), ret_pval = T)
    })
    # print(sprintf("Testing %.2f", v0))
    Pvals = rbind(Pvals, pls)
    t1 = proc.time()[3]
    if(t1 - t0 > 5) {
      print(sprintf("v = %.2f max=%.2f", v0, max(test_vals)))
      t0 = t1
    }
  }

  rownames(Pvals) = round(test_vals, 2)
  return(list(val=test_vals, pvals=Pvals))
}

riTest_hdi = function(y, x, lam0=0) {
  n = length(y); p = ncol(x)

  #LASSO, RIDGE estimates
  beta_lasso = coef(cv.glmnet(x, y, intercept=F))[-1]
  out = solve_ridge(y, x)
  beta_ridge = out$bhat
  e = y - x %*% beta_lasso

  An = out$P_inv %*% t(x)
  tn = function(u) An %*% u

  Tvals = matrix(0, nrow=p)
  print("calculating tvals..")
  for(R in 1:1000) {
    s = sample(c(-1, 1), n, T)
    enew = s * sample(e)  # support heteroskedastic errors.
    Tvals = cbind(Tvals, tn(enew))
  }

  # lam0 = 0
  # Estimate the nuisance term.
  #
  # n0 = out$lam* c(t(lam)%*% out$P_inv %*% beta0)
  # est = tn(e) - bhat + lam0
  # est2 = out$lam* out$P_inv %*% beta_lasso
  # est = out$lam* out$P_inv %*% beta0   ## true value --> great results.
  est =  out$lam* out$P_inv %*% beta_lasso
  # print(sprintf("nuisance=%.2f -- est = %.2f", n0, est))
  print(lam0)

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
  ##
  m = nrow(coef(out))
  # warning("Quick LARS.")
  # return(list(betahat=coef(out)[m/2,]))
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
  # out  = cv.glmnet(x, y, alpha=0, intercept=F)
  warning("Changed lambda definition for ridge regression.")
  p = ncol(x); n = nrow(x)
  lam = sqrt(log(p) /n)

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


main_sim = function(nsim=50, nreps=100, p=500, n=100) {
  cols = c( "pwr_lasso", "pwr_rand", "fwer_lasso", "fwer_rand", "s0", "betas", "HT")
  Results = matrix(0, nrow=0, ncol=length(cols))
  colnames(Results) = cols
  Results = as.data.frame(Results)
  t0 = proc.time()[3]

  Sigma = get_Sigma(p)
  ## s0 = cardinality
  # x = rmvnorm(n, mean=rep(0, p), sigma = Sigma)
  # Z = matrix()
  all_beta_designs = c("U(0,2)", "U(0,4)", "U(-2,2)", "1", "2", "10")

  for(i in 1:nsim) {

    # Sample design. Toeplitz
    x = rmvnorm(n, mean=rep(0, p), sigma = Sigma)
    Z = matrix()
    print(sprintf("Simulation %d/%d...", i, nsim))

    for(s0 in c(3, 15)) {
      # change the support.
      S0 = sample(1:p, size=s0, replace=F)  # true support.
      Snot = setdiff(1:p, S0)

      for(beta_design in all_beta_designs) {
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

        for(ht in c(T, F)) {
          ex = apply(x, 1, function(row) mean(row^2))
          sd_err = (ex * log(1:n+1))^ht; sd_err = sd_err / mean(sd_err)

          ### Now is the main Iteration.
          M = matrix(0, nrow=0, ncol=4)  # interim results.
          colnames(M) = colnames(Results)[1:4]
          ###
          for(irep in 1:nreps) {
            eps = rnorm(n, sd=sd_err) # errors.
            y = x %*% beta + eps

            if(nrow(Z) < 2) {
              lasso = lasso.proj(x, y, return.Z=T)
              Z = lasso$Z
              cache = list(x=x, z=Z)
              save(cache, file=sprintf("out/cache/hdi_cache%d.rda",i))
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
              print(sprintf("iteration %d/%d, s0=%d, beta=%s, HT=%s", irep, nreps, s0, beta_design, ht))
              print(M)
              print(colMeans(M))
              # Store results.
              t0 = t1
            }
          }
          ## Replications over. Save results.
          print("--===-    RESULTS ==--===")
          cm = as.data.frame(t(colMeans(M, na.rm=T)))
          R1 = cbind(cm, data.frame(s0=s0, beta_design=beta_design, ht=ht))
          Results = rbind(Results, R1)
          save(Results, file=sprintf("out/sim_hdi_reps%d_p%d_n%d.rda", nreps, p, n))
          print(Results)
          print("--===- --=-=== ==--===")

        }
      }
    }
  }
}


# Main simulation.
main_sim()