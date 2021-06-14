rm(list=ls())
pacman::p_load(tidyverse)
pacman::p_load(dplyr, stringr)

## Plot summary statistics for RR Tuned and RR Tuning Free
# set working directory
setwd('~/out')
# setwd('~/Dropbox/Research/rr_hdi/code/out/50_100_sign_1')

# set plotting parameters
n <- 50
p <- 100
v <- 'n'
g <- 'sign'

# get csv files in path
files <- list.files(path='.', pattern='.csv', recursive = FALSE)
# separate into CIs and coverages
files_ci <- list.files(path='.', pattern='_ci_quantiles.csv', recursive = FALSE)
files_cov <- setdiff(files, files_ci)

combined_cov <- data.frame()
for (f in 1:length(files_cov)){
  file <- read_csv(files_cov[f])
  file$X1 <- substr(files_cov[f], 1, nchar(files_cov[f]) - 4)
  combined_cov <- rbind(combined_cov, file)
}

combined_ci <- data.frame()
for (f in 1:length(files_ci)){
  file <- read_csv(files_ci[f])
  file$X1 <- substr(files_ci[f], 1, nchar(files_ci[f]) - 17)
  combined_ci <- rbind(combined_ci, file)
}

cov <- paste(v, c("_cov_i", "_cov_a", "_cov_s"), sep = "")
len <- paste(v, c("_len_i", "_len_a", "_len_s"), sep = "")
 
# min value for coverage y-axis
ymin <- .8
# max value for CI length
lenMax <- max(combined_ci$'0.99')

# size of points for coverages
covPch <- -1 # character of points
covCex <- .8 # size of points
covLwd <- 2 # width of line 

# size of points for CIs
maxPch <- 23 # character of point for .99
confIntcex <- .8 # size of point for .99
ciLwd <- 2 # width of line for .25 - .75 quantiles

# x-axis jitter for different methods
rrOffset <- 0.5

# general plotting params
cexAxis_1 <- .6
cexAxis_2 <- .8
sl <- 3.75
sr <- 0.25

# define methods and colors
meth <- c('rr_old_perm_50_100', 
              'rr_perm_50_100')
col <- c('black', 'azure4')

type <- paste(g, n, p, 1, sep = '_')
setEPS()
postscript(paste('../rr_', paste(v, type, sep = '_'), '.eps', sep = ''), width = 7, height = 5)
plot.new()
par(mfcol = c(4,1), oma = c(4.5, 3.0, 1, 0), mar = c(.5, .5, .5, .5))

for (s in c(4, 15)){
  dat <- combined_cov[combined_cov$s == s, ]
  dat <- dat[order(dat$x, dat$e), ]
  dat_q <- combined_ci[combined_ci$s == s, ]
  dat_q <- dat_q[order(dat_q$x, dat_q$e), ]
  
  # deprecated after changes to single_experiment scripts  
  # if ('TGM' %in% dat$x){
  #   dat$x[dat$x == 'TGM'] <- 'GT'
  # }
  # else if ('TG' %in% dat$x){
  #   dat$x[dat$x == 'TG'] <- 'NT'
  # }
  # 
  # if ('HMG' %in% dat$e){
  #   dat$e[dat$e == 'HMG'] <- 'HM'
  # }
  
  labels <- paste(dat$x, dat$e, sep = "-")
  labels_q <- paste(dat_q$x, dat_q$e, sep = "-")
  n_labels <- length(unique(labels))
  
  ## Coverage
  plot((1:length(labels)), dat$a_cov_i, ylim = c(ymin, 1),
       xaxt = "n", ylab = "", xlab = "", type = "n", yaxt = "n",
       xlim = c(sl, (length(labels) - sr)))
  
  # labels for coverage
  if (s == 4){
    mtext(paste("n = ", n, "; p = ", p, "; s = (4, 15)",  sep = ""), cex = .8)
  }
  mtext("Coverage", side = 2, outer = F, line = 2.4, cex = cexAxis_1)
  axis(at = seq(ymin, 1, .05), side = 2, cex.axis = cexAxis_2, las = 2)
  
  # coverage setup
  abline(h = .95, col = "red") ## nominal coverage line
  abline(h = seq(ymin, 1, .1), col = "gray", lty = 3, lwd = .5) ## nominal coverage line
  abline(v = seq(1, length(labels)+2, 1) - rrOffset, col = "darkgray", lty = 2) ## can segment different settings this way
  abline(v = seq(1, length(labels)+2, 0.5) - rrOffset, col = "lightgray", lty = 3) ## can segment different settings this way
  
  # plot coverage
  for (c in 1:length(cov)){
    points((1:n_labels * 2) + 0.5, pull(dat[dat$X1 == meth[1], ][cov[c]]),
           pch = covPch + c, cex = covCex, col = col[1])
    points((1:n_labels * 2) + 0.5 + rrOffset, pull(dat[dat$X1 == meth[2], ][cov[c]]),
           pch = covPch + c, cex = covCex, col = col[2])
  }
  
  ## Quantiles
  plot(-1, -1, xlim = c(sl, (length(labels) - sr)), ylim = c(0, lenMax), xaxt = "n", xlab = "", type = "n", yaxt = "n")
  
  ## labels for quantiles
  mtext("CI Length", side = 2, outer = F, line = 2.4, cex = cexAxis_1)
  axis(side = 2, at = seq(0, lenMax, 2), las = 2, cex.axis = cexAxis_2)
  
  # quantiles setup
  if (s == 15){
    axis(at = c((1:n_labels * 2) - 2 - rrOffset), labels = unique(labels), side = 1, cex.axis = cexAxis_1, las = 2)
  }
  abline(v = seq(1, length(labels)+2, 1) - rrOffset, col = "darkgray", lty = 2) ## can segment different settings this way
  abline(v = seq(1, length(labels)+2, 0.5) - rrOffset, col = "lightgray", lty = 3) ## can segment different settings this way
  abline(h = seq(0, lenMax, 2), col = "lightgray", lty = 3, lwd = .5) ## can segment different settings this way
  
  # plot quantiles
  for (i in 1:length(meth)){
    dat_q_ <- dat_q[(dat_q$X1 == meth[i] & dat_q$Type %in% len), ]
    dat_q_0.25 <- array(0, n_labels)
    dat_q_0.75 <- array(0, n_labels)
    dat_q_0.99 <- array(0, n_labels)
    for (j in 0:(n_labels - 1)){
      dat_q__ <- dat_q_[(j * 3 + 1):((j+1) * 3), ]
      dat_q_0.25[j] <- max(dat_q__$'0.25') 
      dat_q_0.75[j] <- max(dat_q__$'0.75') 
      dat_q_0.99[j] <- max(dat_q__$'0.99') 
    }
    
    segments((1:n_labels * 2) + 0.5 + (i -1) * rrOffset, dat_q_0.25,
             (1:n_labels * 2) + 0.5 + (i -1) * rrOffset, dat_q_0.75,
             lwd = covLwd, col = col[i])
    points((1:n_labels * 2) + 0.5 + (i -1) * rrOffset, dat_q_0.99,
           pch = maxPch, cex = confIntcex, col = col[i])
  }
}

# plot legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("RR Tuned", "RR Tuning Free"), col = c("black", "azure4"),
       pch = 19, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
dev.off()
