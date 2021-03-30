pacman::p_load(hdi, mvtnorm, fastclime, CVXR, Rmosek, torch)

score.nodewiselasso = getFromNamespace("score.nodewiselasso", "hdi")

split_perm_ <- function(n, n_g, X){
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
        XGX_n_m_S[i, , ] <- XGX_n[i, , ] - S
    }
    return(list('XG'=XG, 'XGX_n'=XGX_n, 'XGX_n_m_S'=XGX_n_m_S))
}

min_ell_infinity <- function(XGX_n, S, g)
{
  d <- dim(XGX_n)[1]
  p <- dim(XGX_n)[3]
  # reshape XGX_n_m_S
  XGX_n_ <- aperm(XGX_n, c(1, 3, 2))
  XGX_n <- array(0, c(prod(dim(XGX_n)[1:2]), dim(XGX_n)[3]))
  for (i in 1:d){
    XGX_n[((i - 1) * p + 1):(i * p), ] <- XGX_n_[i, ,]
  }
  
  rm(XGX_n_)
  gc()
  
  # Create constraints for t
  XGX_constraints <- array(0, c(prod(dim(XGX_n)[1]), d+1))
  for (i in 1:d){
    XGX_temp <- array(0, c(p, d+1))
    XGX_temp[, i] <- 1
    XGX_constraints[((i - 1) * p + 1):(i * p), ] <- XGX_temp
  }
  S_constraints <- array(0, c(p, d+1))
  S_constraints[, d+1] <- 1
  # S_constraints <- array(0, c(p, 1))
  # S_constraints[, 1] <- 1
  # Linear constraints in [x; t]
  prob <- list(sense="min")
  prob$A <- rbind(cbind(XGX_n, -XGX_constraints),
                  cbind(XGX_n, XGX_constraints),
                  cbind(-S, -S_constraints),
                  cbind(-S, S_constraints))
  # prob$A <- rbind(cbind(-S, -S_constraints),
  #                 cbind(-S, S_constraints))
  # Bound values for constraints
  prob$bc <- rbind(blc=c(array(rep(-Inf, dim(XGX_n)[1])), array(rep(0, dim(XGX_n)[1])), rep(-Inf, p), -g),
                   buc=c(array(rep(0, dim(XGX_n)[1])), array(rep(Inf, dim(XGX_n)[1])), -g, rep(Inf, p)))
  # prob$bc <- rbind(blc=c(rep(-Inf, p), -g),
  #                  buc=c(-g, rep(Inf, p)))
  # Bound values for variables
  prob$bx <- rbind(blx=c(rep(-Inf, p), rep(0, d+1)),
                   bux=rep(Inf, p+d+1))
  # prob$bx <- rbind(blx=c(rep(-Inf, p), rep(0, 1)),
  #                  bux=rep(Inf, p+1))
  # Coefficients for variables
  prob$c <- c(rep(0, p), 0.01 * rep(1 / d,  d), 0.99)
  # prob$c <- c(rep(0, p), 1)
  # prob$sol <- list(bas=list())
  # prob$sol$bas$xx <- c(m_i, rep(0, d+1))
  
  # prob$iparam <- list(OPTIMIZER="OPTIMIZER_INTPNT")
  prob$iparam <- list(OPTIMIZER="OPTIMIZER_FREE_SIMPLEX")
  r <- mosek(prob, list(verbose=1))
  # return(list('m'=r$sol$itr$xx[1:p]))
  m <- r$sol$bas$xx[1:p]
  # print(r$sol$bas$xx[p+1])
  rm(prob, r)
  gc()
  
  return(list('m'=m))
  
}


set.seed(1)

n <- 50
p <- 100
n_g <- 500
n_e <- 50
idx <- 1

fc <- rep(0, n_e)
to <- rep(0, n_e)
# l2 <- rep(0, n_e)
# gu <- rep(0, n_e)
cvx_msk <- rep(0, n_e)
# scs <- rep(0, n_e)
diff <- rep(0, n_e)
time <- rep(0, n_e)


for (i in 1:n_e){
    x <- mvtnorm::rmvnorm(n, mean=rep(0, p), sigma=diag(p))
    S <- t(x) %*% (x) / n

    group_actions <- split_perm_(n, n_g, x)
    XGX_n_m_S <- group_actions$XGX_n_m_S
    XGX_n <- group_actions$XGX_n
    # 
    # XGX_n <- torch_tensor(group_actions$XGX_n)
    I <- diag(p)

    i_t <- rep(0, p)
    i_t[idx] <- 1
    # 
    # # fastclime
    # lambda <- 0.1 * sqrt(log(p) / n)
    # # Get path of Ms from fastclime
    # clime_M <- fastclime(S, lambda.min=lambda, nlambda=500)
    # M <- clime_M$icovlist[[length(clime_M$icovlist)]]
    # m_ <- M[, idx]
    # 
    # tau <- 0
    # for (j in 1:n_g){
    #     tau <- tau + max(abs(t(XGX_n_m_S[j, , ]) %*% m_ + i_t ))
    # }
    # fc[i] <- tau / n_g
    # print(mean(fc[1:i]))
    # rm(clime_M)
    # gc()

    start <- proc.time()
    msk <- min_ell_infinity(XGX_n, S, i_t)
    tau <- 0
    for (j in 1:n_g){
      tau <- tau + max(abs(t(XGX_n[j, , ]) %*% msk$m + i_t ))
    }
    # print(tau / n_g)
    time[i] <- proc.time() - start
    
    # # nodewise
    # node <- score.nodewiselasso(x, wantTheta = TRUE, verbose = FALSE,
    #                             lambdaseq = "quantile", parallel = FALSE, ncores = 4,
    #                             oldschool = FALSE, lambdatuningfactor = 1)
    # M <- node$out
    # m <- M[, idx]
    #
    # tau <- 0
    # for (j in 1:n_g){
    #     tau <- tau + max(abs(t(XGX[[j]] - S) %*% m + i_t ))
    # }
    # nl[i] <- tau / n_g
    # print(mean(nl[1:i]))
    # rm(node)
    # gc()
    
    # # Get M
    # # Setup
    # tau <- Variable(2)
    # beta <- 50
    
    # beta <- 200
    # learning_rate <- 1e-3
    # iterations <- 10000
    # 
    # S_ <- torch_tensor(S)
    # 
    # loss_approx <- function(idx, m_i, beta) {
    #     # return(torch_logsumexp(beta * torch_abs(I[, idx] - torch_matmul(S_, m_i)), dim = 1) / beta +
    #     #            torch_mean(torch_logsumexp(beta * torch_abs(torch_matmul(XGX_n, m_i)), dim = 2) / beta))
    #     return(
    #         torch_mean(
    #         torch_logsumexp(
    #             beta * torch_abs(I[, idx] + torch_matmul(torch_transpose(XGX_n - S_, 3, 2), m_i)), dim = 2
    #             ) / beta)
    #         )
    # }
    # 
    # for (idx in 1){
    #     # Initialise variables
    #     # m_i <- torch_zeros(p, requires_grad = TRUE)
    #     m_i <- torch_tensor(m_, requires_grad = TRUE)
    #     
    #     for (iter in 1:iterations){
    #         loss <- loss_approx(1, m_i, beta)
    #         loss$backward()
    #         
    #         with_no_grad({
    #             m_i <- m_i$sub_(learning_rate * m_i$grad)
    #             m_i$grad$zero_()
    #         })
    #     }
    # }
    # 
    # tau <- 0
    # for (j in 1:n_g){
    #     tau <- tau + max(abs(t(XGX_n_m_S[j, , ]) %*% as_array(m_i) + i_t ))
    # }
    # m_to <- as_array(m_i)
    # to[i] <- tau / n_g
    # print(mean(to[1:i]))
    
    
    # start <- proc.time()
    # 
    # XGX_n_m_S <- group_actions$XGX_n_m_S
    # # XGX_n <- group_actions$XGX_n
    # I_i <- rep(0, p)
    # I_i[idx] <- 1
    # # tau_1 <- Variable(1)
    # tau_2 <- Variable(n_g)
    # tau_2@value <- rep(0, n_g)
    # m_i <- Variable(p)
    # m_i@value <- rep(0, p)
    # # m_i@value <- m_
    # # 
    # obj <- Minimize(mean(tau_2))
    # # obj <- Minimize(mean(tau_1))
    # # obj <- Minimize(tau_1 + 1 / n_g * sum_entries(tau_2))
    # # # cost <- 0
    # # # for (j in 1:n_g){
    # # #     cost <- cost + t(XGX[[j]] - S) %*% m + i_t
    # # #     # constraints <- c(constraints, t(XGX[[j]] - S) %*% m + i_t >= -tau * rep(1, p))
    # # #     # constraints <- c(constraints, XGX[[j]] %*% t(M_[i, ]) <= lambda * matrix(1, p, 1))
    # # #     # constraints <- c(constraints, XGX[[j]] %*% t(M_[i, ]) >= -lambda * matrix(1, p, 1))
    # # # }
    # # # obj <- Minimize(cvxr_norm(cost))
    # # # prob <- Problem(obj)
    # # # start <- proc.time()
    # # # result <- solve(prob, solver="GUROBI")
    # # # tau <- 0
    # # # for (j in 1:n_g){
    # # #     tau <- tau + max(abs(t(XGX[[j]] - S) %*% result$getValue(m) + i_t ))
    # # # }
    # # # l2[i] <- tau / n_g
    # # # print(proc.time() - start)
    # # # print(mean(l2[1:i]))
    # # # rm(result)
    # # # gc()
    # # 
    # # # browser()
    # # 
    # # Define constraints
    # constraints <- list()
    # for (j in 1:n_g){
    #   # constraints <- c(constraints, t(XGX_n[j, , ]) %*% m_i <= tau_2[j] * rep(1, p))
    #   # constraints <- c(constraints, t(XGX_n[j, , ]) %*% m_i >= -tau_2[j] * rep(1, p))
    #   constraints <- c(constraints, t(XGX_n_m_S[j, , ]) %*% m_i + I_i <= tau_2[j] * rep(1, p))
    #   constraints <- c(constraints, t(XGX_n_m_S[j, , ]) %*% m_i + I_i >= -tau_2[j] * rep(1, p))
    # }
    # # constraints <- c(constraints, - S %*% m_i + I_i <= tau_1 * rep(1, p))
    # # constraints <- c(constraints, - S %*% m_i + I_i >= -tau_1 * rep(1, p))
    # 
    # prob <- Problem(obj, constraints)
    # # # gurobi
    # # start <- proc.time()
    # # tau <- 0
    # # result <- solve(prob, solver="GUROBI")
    # # for (j in 1:n_g){
    # #     tau <- tau + max(abs(t(XGX_n_m_S[j, , ]) %*% result$getValue(m) + i_t ))
    # # }
    # # rm(result)
    # # gc()
    # # 
    # # # gu_solve <- function(){
    # # #     result <- solve(prob, solver="GUROBI")
    # # #     for (j in 1:n_g){
    # # #         tau <- tau + max(abs(t(XGX[[j]] - S) %*% result$getValue(m) + i_t ))
    # # #     }
    # # #     rm(result)
    # # #     gc()
    # # # }
    # # # 
    # # # try(gu_solve(), silent = TRUE)
    # # 
    # # gu[i] <- tau / n_g
    # # print(mean(gu[1:i]))
    # # print(proc.time() - start)
    # # 
    # # mosek
    # # start <- proc.time()
    # solver_opts <- list()
    # # solver_opts$iparam <- list(OPTIMIZER = "MSK_OPTIMIZER_FREE_SIMPLEX")
    #     # list(MSK_IPAR_OPTIMIZER="MSK_OPTIMIZER_FREE_SIMPLEX")
    # tau <- 0
    # result <- solve(prob, solver="MOSEK", solver_opts=solver_opts)
    # for (j in 1:n_g){
    #     tau <- tau + max(abs(t(XGX_n_m_S[j, , ]) %*% result$getValue(m_i) + i_t ))
    # }
    # m_msk <- result$getValue(m_i) 
    # rm(result)
    # gc()
    # 
    # # msk_solve <- function(){
    # #     result <- solve(prob, solver="MOSEK")
    # #     for (j in 1:n_g){
    # #         tau <- tau + max(abs(t(XGX[[j]] - S) %*% result$getValue(m) + i_t ))
    # #     }
    # #     rm(result)
    # #     gc()
    # # }
    # # 
    # # try(msk_solve(), silent = TRUE)
    # 
    # cvx_msk[i] <- tau / n_g
    # print(mean(cvx_msk[1:i]))
    # print(proc.time() - start)
    # 
    # # # scs
    # # start <- proc.time()
    # # tau <- 0
    # # result <- solve(prob, solver="SCS")
    # # for (j in 1:n_g){
    # #     tau <- tau + max(abs(t(XGX_n_m_S[j, , ]) %*% result$getValue(m) + i_t ))
    # # }
    # # rm(result)
    # # gc()
    # # 
    # # # scs_solve <- function(){
    # # #     result <- solve(prob, solver="SCS")
    # # #     for (j in 1:n_g){
    # # #         tau <- tau + max(abs(t(XGX[[j]] - S) %*% result$getValue(m) + i_t ))
    # # #     }
    # # #     rm(result)
    # # #     gc()
    # # # }
    # # # 
    # # # try(scs_solve(), silent = TRUE)
    # # 
    # # scs[i] <- tau / n_g
    # # print(mean(scs[1:i]))
    # # print(proc.time() - start)
    # # 
    # rm(prob, constraints)
    # diff[i] <- sqrt(sum((msk$m - m_msk)**2))
    # browser()
}

print(fc)
# print(mean(nl))
print(gu)
print(msk)
print(scs)

# max(abs(diag(p) - node$out %*% t(x) %*% x / n))

# sd <- attr(scale(x), "scaled:scale")
# me <- attr(scale(x), "scaled:center")
# x_ <- scale(x, center = FALSE)
# S <- t(x_) %*% (x_) / n
# S_bar <- t(x_ - colMeans(x_)) %*% (x_ - colMeans(x_)) / n
# # Solve for M
# lambda <- 0.01 * sqrt(log(p) / n)
# # Get path of Ms from fastclime
# clime_M <- fastclime(S, lambda.min=lambda, nlambda=500)
#
# max(abs(diag(p) - clime_M$icovlist[[length(clime_M$icovlist)]] %*% S))
# norm(clime_M$icovlist[[length(clime_M$icovlist)]], '1')
#
# beta <- sample(c(-1, 1), size=p, replace=T)
# err <- rmvnorm(1, mean=rep(0, n), sigma=diag(n))
#
# y = x %*% beta + as.vector(err)
#
# x_ %*% ((beta * sd))
# x %*% beta