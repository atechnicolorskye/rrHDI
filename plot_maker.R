pacman::p_load(dplyr, stringr)

dat_raw <- read.csv("out/combined_means_hard.csv")
dat_quantiles <- read.csv("out/combined_quantiles_hard.csv")

n <- 50
p <- 100
v <- 'n'
g <- 'perm'
cov <- paste(v, 'cov', sep = '_')
type <- paste(g, n, p, sep = '_')

# min value for coverage y-axis
# ymin <- .0
ymin <- .5
# max value for CI length
lenMax <- max(dat_quantiles[str_detect(dat_quantiles$Type, type), 'X0.99'])
# lenMax <- 6.5
# size of points for CI's

medianPch <- "" # character of median
medianCex <- 1.8 # size of point for median
maxPch <- 19 # character of point for .99
confIntcex <- .3 # size of point for .99
ciLwd <- 2.5 # width of line for .25 - .75 quantiles

## x-axis jitter for different methods
rrOffset <- .3
hdiOffset <- .15
blprOffset <- 0
dlassoOffset <- -.15
silmOffset <- -.3

setEPS()
postscript(paste(paste(v, type, sep = '_'), '_hard.eps', sep = ''), width = 7, height = 4)
# postscript(paste(paste(v, type, sep = '_'), '.eps', sep = ''), width = 7, height = 4)
par(mfcol = c(2,2), oma = c(4.5, 2.5, 1, 0), mar = c(.5, .5, .5, .5))

for (s in c(3, 10)){
  dat <- dat_raw[str_detect(dat_raw$Type, type), ]
  dat <- dat[dat$s == s, ]
  dat <- dat[order(dat$x, dat$e), ]
  dat_q <- dat_quantiles[str_detect(dat_quantiles$Type, type), ]
  dat_q <- dat_q[dat_q$s == s, ]
  dat_q <- dat_q[order(dat_q$x, dat_q$e), ]
  labels <- paste(dat$x, dat$e, sep = "-")
  labels_q <- paste(dat_q$x, dat_q$e, sep = "-")
  n_labels <- length(unique(labels))
  
  # Coverage
  plot((1:length(labels) - rrOffset), select(dat, cov)[, 1], ylim = c(ymin, 1),
       xaxt = "n", ylab = "", xlab = "", pch = 19, cex = .6, yaxt = "n",
       xlim = c(1, length(labels)))
  
  ## Labels for coverage
  if (s == 3){
    mtext("Coverage", side = 2, outer = F, line = 2.2, cex = .7)
    axis(at = seq(ymin, 1, .1), side = 2, cex.axis = .8, las = 2)
  } 
  mtext(paste("n = ", n, "; p = ", p, "; s = ",s, sep = ""))
  
  # Coverage setup
  abline(h = .95, col = "red") ## nominal coverage line
  abline(h = seq(ymin, 1, .1), col = "gray", lty = 3, lwd = .5) ## nominal coverage line
  abline(v = seq(1, length(labels)+1, 5) - rrOffset, col = "darkgray", lty = 2) ## can segment different settings this way
  abline(v = seq(1, length(labels), 1) - rrOffset, col = "lightgray", lty = 3) ## can segment different settings this way
  
  # Coverage proper
  points(((1:n_labels * 5) - 4 - rrOffset), select(dat[str_detect(dat$Type, 'rr'), ], cov)[, 1], col = "black", pch = 19, cex = .6)
  points(((1:n_labels * 5) - 3 - rrOffset), select(dat[str_detect(dat$Type, 'hdi'), ], cov)[, 1], col = "blue3", pch = 19, cex = .6)
  points(((1:n_labels * 5) - 2 - rrOffset), select(dat[str_detect(dat$Type, 'blpr'), ], cov)[, 1], col = "chartreuse4", pch = 19, cex = .6)
  points(((1:n_labels * 5) - 1 - rrOffset), select(dat[str_detect(dat$Type, 'dlasso'), ], cov)[, 1], col = "darkorange", pch = 19, cex = .6)
  points(((1:n_labels * 5) - 0 - rrOffset), select(dat[str_detect(dat$Type, 'silm'), ], cov)[, 1], col = "deeppink3", pch = 19, cex = .6)
  
  ## Quantiles
  plot(-1, -1, xlim = c(1, length(labels)), ylim = c(0, lenMax), xaxt = "n", xlab = "", pch = 19, cex = .6, yaxt = "n")
  
  ## Labels for quantiles
  if (s == 3){
    axis(side = 2, at = seq(0, lenMax, 2), las = 2, cex.axis = .8)
  mtext("CI Length", side = 2, outer = F, line = 2.2, cex = .7)
  }
  
  ## Quantiles setup
  axis(at = c((1:n_labels * 5) - 2 - rrOffset), labels = unique(labels), side = 1, cex.axis = .8, las = 2)
  abline(v = seq(1, length(labels)+1, 5) - rrOffset, col = "darkgray", lty = 2) ## can segment different settings this way
  abline(v = seq(1, length(labels), 1) - rrOffset, col = "lightgray", lty = 3) ## can segment different settings this way
  abline(h = seq(0, lenMax, 2), col = "lightgray", lty = 3, lwd = .5) ## can segment different settings this way
  
  ## Quantiles proper
  ## Small
  # segments((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_4000\\b', sep = '_')), "X0.25"],
  #          (1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_4000\\b', sep = '_')), "X0.75"],
  #          lwd = ciLwd)
  # points((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_4000\\b', sep = '_')), c("X0.5")],
  #        pch = medianPch, cex = medianCex)
  # points((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_4000\\b', sep = '_')), c("X0.99")],
  #        pch = maxPch, cex = confIntcex)
  # # Medium
  segments((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_10000\\b', sep = '_')), "X0.25"],
           (1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_10000\\b', sep = '_')), "X0.75"],
           lwd = ciLwd)
  points((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_10000\\b', sep = '_')), c("X0.5")],
         pch = medianPch, cex = medianCex)
  points((1:n_labels * 5) - 4 - rrOffset, dat_q[str_detect(dat_q$Type, paste('\\brr', type, v, 'len_rr_10000\\b', sep = '_')), c("X0.99")],
         pch = maxPch, cex = confIntcex)

  segments((1:n_labels * 5) - 3 - rrOffset, dat_q[str_detect(dat_q$Type, paste('hdi', type, v, 'len', sep = '_')), "X0.25"],
           (1:n_labels * 5) - 3 - rrOffset, dat_q[str_detect(dat_q$Type, paste('hdi', type, v, 'len', sep = '_')), "X0.75"], 
           lwd = ciLwd, col = "blue3")
  points((1:n_labels * 5) - 3 - rrOffset, dat_q[str_detect(dat_q$Type, paste('hdi', type, v, 'len', sep = '_')), c("X0.5")], 
         pch = medianPch, cex = medianCex, col = "blue3")
  points((1:n_labels * 5) - 3 -rrOffset, dat_q[str_detect(dat_q$Type,paste('hdi', type, v, 'len', sep = '_')), c("X0.99")], 
         pch = maxPch, cex = confIntcex, col = "blue3")
    
  segments((1:n_labels * 5) - 2 - rrOffset, dat_q[str_detect(dat_q$Type, paste('blpr', type, v, 'len', sep = '_')), "X0.25"],
           (1:n_labels * 5) - 2 - rrOffset, dat_q[str_detect(dat_q$Type, paste('blpr', type, v, 'len', sep = '_')), "X0.75"], 
           lwd = ciLwd, col = "chartreuse4")
  points((1:n_labels * 5) - 2 - rrOffset, dat_q[str_detect(dat_q$Type, paste('blpr', type, v, 'len', sep = '_')), c("X0.5")], 
         pch = medianPch, cex = medianCex, col = "chartreuse4")
  points((1:n_labels * 5) - 2 - rrOffset, dat_q[str_detect(dat_q$Type, paste('blpr', type, v, 'len', sep = '_')), c("X0.99")], 
         pch = maxPch, cex = confIntcex, col = "chartreuse4")
  
  segments((1:n_labels * 5) - 1 - rrOffset, dat_q[str_detect(dat_q$Type, paste('dlasso', type, v, 'len', sep = '_')), "X0.25"],
           (1:n_labels * 5) - 1 - rrOffset, dat_q[str_detect(dat_q$Type, paste('dlasso', type, v, 'len', sep = '_')), "X0.75"], 
           lwd = ciLwd, col = "darkorange")
  points((1:n_labels * 5) - 1 - rrOffset, dat_q[str_detect(dat_q$Type, paste('dlasso', type, v, 'len', sep = '_')), c("X0.5")], 
         pch = medianPch, cex = medianCex, col = "darkorange")
  points((1:n_labels * 5) - 1 - rrOffset, dat_q[str_detect(dat_q$Type, paste('dlasso', type, v, 'len', sep = '_')), c("X0.99")], 
         pch = maxPch, cex = confIntcex, col = "darkorange")
  
  segments((1:n_labels * 5) - 0 - rrOffset, dat_q[str_detect(dat_q$Type, paste('silm', type, v, 'len', sep = '_')), "X0.25"],
           (1:n_labels * 5) - 0 - rrOffset, dat_q[str_detect(dat_q$Type, paste('silm', type, v, 'len', sep = '_')), "X0.75"], 
           lwd = ciLwd, col = "deeppink3")
  points((1:n_labels * 5) - 0 - rrOffset, dat_q[str_detect(dat_q$Type, paste('silm', type, v, 'len', sep = '_')), c("X0.5")], 
         pch = medianPch, cex = medianCex, col = "deeppink3")
  points((1:n_labels * 5) - 0 - rrOffset, dat_q[str_detect(dat_q$Type, paste('silm', type, v, 'len', sep = '_')), c("X0.99")], 
         pch = maxPch, cex = confIntcex, col = "deeppink3")
}

# Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("RR", "HDI", "BLPR", "DLASSO", "SILM"), col = c("black","blue3", "chartreuse4", "darkorange", "deeppink3"),
       pch = 19, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')

dev.off()
