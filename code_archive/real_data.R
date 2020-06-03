rm(list=ls())
# Libraries
library("doParallel")
library("fastclime")
library("glmnet")
library("mvtnorm")
library("HDCI") # Buhlmann
# Code
source("randomization_hdi.R") # Us
# Data
library("care")

seed <- 2020
set.seed(seed)

data(efron2004)
n <- dim(efron2004$x)[0]
p <- dim(efron2004$x)[1]

x <- matrix(efron2004$x, c(n, p))
y <- matrix(efron2004$y)

ind <- 1:10
beta_0 <- rep(0, 10)

S <- (t(x) %*% x) / n
lambda_min <- min(.1, sqrt(log(p) / n))
n_solve <- 50
dantzig_M <- fastclime(x, lambda.min=lambda_min, nlambda=n_solve)
out = riSingle_hdi_dantzig(y, x, dantzig_M$sigmahat, S_test=ind, lam0_vals=beta_0)
