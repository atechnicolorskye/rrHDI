pacman::p_load(dplyr, stringr)

dat_raw <- read.csv("out/combined_means_silm_rr.csv")

n <- 50
p <- 100
v <- 'a'
g <- 'sign'
cov <- paste(v, 'cov', sep = '_')
type <- paste(g, n, p, sep = '_')

# min value for coverage y-axis
ymin <- .8
# ymin <- .5
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
rrOffset <- .05
silmOffset <- -.05

setEPS()
postscript('silm_rr_small_sign.eps', width = 7, height = 3)
# postscript(paste(paste(v, type, sep = '_'), '.eps', sep = ''), width = 7, height = 4)
par(mfcol = c(1,2), oma = c(4.5, 2.5, 1, 0), mar = c(.5, .25, .5, .25))

for (s in c(3, 10)){
    dat <- dat_raw[str_detect(dat_raw$Type, type), ]
    dat <- dat[dat$s == s, ]
    dat <- dat[order(dat$x, dat$e), ]
    labels <- paste(dat$x, dat$e, sep = "-")
    n_labels <- length(unique(labels))
    
    # Coverage
    plot(-1, -1, ylim = c(ymin, 1),
         xaxt = "n", ylab = "", xlab = "", pch = 19, cex = .6, yaxt = "n",
         xlim = c(1.5, length(labels)))
    
    ## Labels for coverage
    if (s == 3){
        mtext("Coverage", side = 2, outer = F, line = 2, cex = .7)
        axis(at = seq(ymin, 1, .05), side = 2, cex.axis = .6, las = 2)
    }
    mtext(paste("n = ", n, "; p = ", p, "; s = ",s, sep = ""))
    
    # Coverage setup
    abline(h = .95, col = "red") ## nominal coverage line
    abline(h = seq(ymin, 1, .1), col = "gray", lty = 3, lwd = .5) ## nominal coverage line
    abline(v = seq(1, length(labels)+1, 2) - rrOffset, col = "darkgray", lty = 2) ## can segment different settings this way
    abline(v = seq(1, length(labels), 1) - rrOffset, col = "lightgray", lty = 3) ## can segment different settings this way
    
    # Coverage proper
    points(((1:n_labels * 2) - 1 - rrOffset), select(dat[str_detect(dat$Type, 'rr'), ], cov)[, 1], col = "black", pch = 19, cex = .6)
    # points(((1:n_labels * 5) - 3 - rrOffset), select(dat[str_detect(dat$Type, 'hdi'), ], cov)[, 1], col = "blue3", pch = 19, cex = .6)
    # points(((1:n_labels * 5) - 2 - rrOffset), select(dat[str_detect(dat$Type, 'blpr'), ], cov)[, 1], col = "chartreuse4", pch = 19, cex = .6)
    # points(((1:n_labels * 5) - 1 - rrOffset), select(dat[str_detect(dat$Type, 'dlasso'), ], cov)[, 1], col = "darkorange", pch = 19, cex = .6)
    points(((1:n_labels * 2) - 0 - rrOffset), select(dat[str_detect(dat$Type, 'silm'), ], cov)[, 1], col = "deeppink3", pch = 19, cex = .6)
    
    ## Quantiles setup
    axis(at = c((1:n_labels * 2) - rrOffset), labels = unique(labels), side = 1, cex.axis = .7, las = 2)
}

# Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("RR", "SILM"), col = c("black", "deeppink3"),
       pch = 19, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')

dev.off()
