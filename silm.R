library(SILM) # Zhang and Guang

# Fixed SILM code
SimE.CI = function (X, Y, a_ind, n_ind, M = 1000, alpha = 0.95) 
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
  
  # Active set
  a_band.nst <- matrix(NA, length(a_ind), 2)
  for (i in 1:length(a_ind)){
    a_stat.boot.nst <- rep(NA, M)
    for (j in 1:M) {
      e <- rnorm(n)
      xi.boot <- Theta[a_ind[i], ] %*% t(X) %*% e * sqrt(sigma.sq)/sqrt(n)
      a_stat.boot.nst[j] <- max(abs(xi.boot))
    }
    a_crit.nst <- quantile(a_stat.boot.nst, alpha)
    a_band.nst[i, ] <- cbind(beta.db[a_ind[i]] - a_crit.nst/sqrt(n), beta.db[a_ind[i]] + a_crit.nst/sqrt(n))
  }
  
  # Non-active set
  n_band.nst <- matrix(NA, length(n_ind), 2)
  for (i in 1:length(n_ind)){
    n_stat.boot.nst <- rep(NA, M)
    for (j in 1:M) {
      e <- rnorm(n)
      xi.boot <- Theta[n_ind[i], ] %*% t(X) %*% e * sqrt(sigma.sq)/sqrt(n)
      n_stat.boot.nst[j] <- max(abs(xi.boot))
    }
    n_crit.nst <- quantile(n_stat.boot.nst, alpha)
    n_band.nst[i, ] <- cbind(beta.db[n_ind[i]] - n_crit.nst/sqrt(n), beta.db[n_ind[i]] + n_crit.nst/sqrt(n))
  }
  
  result <- list(beta.db[c(a_ind, n_ind)], a_band.nst, n_band.nst)
  names(result) <- c("de-biased Lasso", "a_band.nst", "n_band.nst")
  return(result)
}