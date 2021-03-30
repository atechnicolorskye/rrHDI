split_perm <- function(n, n_g, X){
  XG <- array(0, dim=c(n_g, dim(X)[2], dim(X)[1]))
  XGX_n <- array(0, dim=c(n_g, dim(X)[2], dim(X)[2]))
  XGX_n_m_S <- array(0, dim=c(n_g, dim(X)[2], dim(X)[2]))
  
  S <- t(X) %*% X / n
  
  fir <- sample(n, n/2)
  sec <- setdiff(1:n, fir)
  for (i in 1:n_g){
    ind <- rep(0, n)
    ind[fir] <- sample(sec)
    ind[sec] <- sample(fir)
    G <- diag(length(ind))[ind,]
    XG[i, , ] <- t(X) %*% G
    XGX_n[i, , ] <- XG[i, , ] %*% X / n
  }
  
  return(list('XG'=XG, 'XGX_n'=XGX_n, 'XGX_n_m_S'=XGX_n_m_S))
}

### Function to calculate
d_lambda <- function(M, xgx, S, a){
  p <- ncol(M)
  ### calculates the two terms broken apart ### 
  t1 <- diag(p)[a, , drop = F] - t(M)[a, , drop = F ] %*% S
  t2 <- apply(xgx$XGX_n, MAR = 1, function(z){t(M)[a, , drop = F] %*% z})
  d1 <- max(abs(t1)) + mean(apply(t2, MAR = 2, function(z){max(abs(z))}))
                                   
                                  
  ## calculates two terms together
  d2 <- mean(unlist(apply(t2, MAR = 2, function(z){max(abs(t1 + z))})))
  
  return(list(d1 = d1, d2 = d2, bias = max(abs(t1))))
}



### Set up problem ###
n <- 100
p <- 300
M <- toeplitz(.8^c(0:(p-1)))
M <- diag(p)
X <- mvtnorm::rmvnorm(n, sigma = M)
xgx <- split_perm(n, 200, X)
S <- t(X) %*% X/n

a <- 1

l1 <- 2 * sqrt(log(p)/n)
system.time(climeOutOld <- clime::clime(S, sigma = T, lambda = seq(1, l1, by = -.2)))

MInit <- climeOutOld$Omegalist[[length(climeOutOld$Omegalist)]]
max(abs(S %*% MInit - diag(p)))

climeOutInit <- fastclime::fastclime(S, lambda.min =  .1, nlambda =50)
icovListInit <- sapply(seq(1, .1, by = -.05),
                   FUN = function(lam){fastclime::fastclime.selector(climeOutInit$lambdamtx, climeOutInit$icovlist, lambda = lam)$icov},
                   simplify = F)

dInit <- d_lambda(MInit, xgx, S, a)

climeOut <- fastclime::fastclime(S, lambda.min = .1, nlambda =150)
icovList <- sapply(climeOut$lambdamtx[, a],
                   FUN = function(lam){fastclime::fastclime.selector(climeOut$lambdamtx, climeOut$icovlist, lambda = lam)$icov},
                   simplify = F)
### theoretically driven M ###
d1 <- d_lambda(MInit, xgx, S, a)
d_list <- sapply(icovList, d_lambda, xgx, S, a)


par(mfrow = c(1,3))
plot(unlist(d_list[1, ]))
abline(h = d1[1])
plot(unlist(d_list[2, ]))
abline(h = d1[2])
plot(unlist(d_list[3, ]))
abline(h = d1[3])


