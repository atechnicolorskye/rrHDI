rm(list=ls())
library(tidyverse)

## Compile results into single csv file
# setwd('out')
# filenames <- list.files(pattern = "_n_50_", full.names = TRUE)
# df <- filenames %>% lapply(read_csv) %>%  bind_rows
# df <- df[, -1]
# write.csv(df, 'small.csv')

# filenames <- list.files(pattern = "n_100_", full.names = TRUE)
# df <- filenames %>% lapply(read_csv) %>%  bind_rows
# df <- df[, -1]
# write.csv(df, 'large.csv')


library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

## Process results
# *_admissible are results filtered for admissible group actions
# *_admissible_drop are results without unscaled RR
data <- read.csv("out/b_1000_n_solve_500/large_admissible_drop.csv", stringsAsFactors = F)
data[, 1:20] <- as.numeric(as.matrix(data[, 1:20]))

win.rateA <- tabulate(apply(data[, 1:5], MAR = 1, function(x){which.min(abs(x - .95))}), nbins=5)
win.rateN <- tabulate(apply(data[, 11:15], MAR = 1, function(x){which.min(abs(x - .95))}), nbins=5)

# Min, Q2, Q3, Max
tab1 <- round(cbind(apply(data[, 1:5], MAR = 2, min, na.rm =T),
                    apply(data[, 1:5], MAR = 2, median, na.rm =T),
                    apply(data[, 1:5], MAR = 2, function(x){quantile(x, .75)}),
                    apply(data[, 1:5], MAR = 2, max, na.rm =T),
                    apply(data[, 11:15], MAR = 2, min, na.rm =T),
                    apply(data[, 11:15], MAR = 2, median, na.rm =T),
                    apply(data[, 11:15], MAR = 2, function(x){quantile(x, .75)}),
                    apply(data[, 11:15], MAR = 2, max, na.rm =T)), 3)

colnames(tab1) <- c("min CR-A", "Q2 CR-A", "Q3 CR-A", "Max CR-A", "min CR-N", "Q2 CR-N", "Q3 CR-N", "Max CR-N")
rownames(tab1) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")

# Various statistics
tab2 <- round(cbind(apply(abs(data[, 1:5] - .95), MAR = 2, mean, na.rm =T),
                    apply(abs(data[, 1:5] - .95), MAR = 2, median, na.rm =T),
                    apply(data[, 6:10], MAR = 2, mean, na.rm =T),
                    apply(data[, 6:10], MAR = 2, median, na.rm =T),
                    win.rateA,
                    apply(abs(data[, 11:15] - .95), MAR = 2, mean, na.rm =T),
                    apply(abs(data[, 11:15] - .95), MAR = 2, median, na.rm =T),
                    apply(data[, 16:20], MAR = 2, mean, na.rm =T),
                    apply(data[, 16:20], MAR = 2, median, na.rm =T),
                    win.rateN), 3)

colnames(tab2) <- c("Mn CR-A", "Md CR-A", "Mn LN-A", "Md LN-A", "Win-A", "Mn CR-N", "Md CR-N", "Mn LN-N", "Md LN-N", "Win-N")
rownames(tab1) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")

# Rankings
rank_dist_a <- apply(data[, 1:5], MAR = 1, function(x){abs(x - .95)})
tab3 <- t(apply(rank_dist_a, MAR = 2, rank))
colnames(tab3) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")
rownames(tab3) <- paste(data$s, data$x, data$b, data$e, data$g, sep=', ')
means <- round(apply(tab3, MAR = 2, mean), 2)
tab3 <- rbind(tab3, means)
rownames(tab3)[rownames(tab3) == "means"] <- "Average Rank for Active Parameters"

xtable(tab3)

rank_dist_n <- apply(data[, 11:15], MAR = 1, function(x){abs(x - .95)})
tab4 <- t(apply(rank_dist_n, MAR = 2, rank))
colnames(tab4) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")
rownames(tab4) <- paste(data$s, data$x, data$b, data$e, data$g, sep=', ')
means <- round(apply(tab4, MAR = 2, mean), 2)
tab4 <- rbind(tab4, means)
rownames(tab4)[rownames(tab4) == "means"] <- "Average Rank for Inactive Parameters"

xtable(tab4)

reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

# Plots
pdf('small.pdf')
num.sim.settings <- dim(data)[1]
sim.per.row <- 10
rows <- ceiling(num.sim.settings / sim.per.row)

par(mfrow = c(rows, sim.per.row), mar = c(0, .1, 1.5, 0), oma = c(3, 2.1, 0, .3))
ind <- 1
for(j in 1:rows){
  i <- 1
  
  plot.char <- rep(19, 5)
  plot.char[which.min(abs(data[ind, 1:5] - .95))] <- 1  
  plot(-1, -1, xlim = c(0, 1), ylim = c(.7, 1), xlab = "" ,ylab = "",
       xaxt = "n")
  mtext(text = paste((data[ind, 21:25]), collapse = "."), line = 0, cex = .5)
  abline(h = seq(.5, 1, by = .05), lty = 2, col = "gray", lwd = .5)
  abline(h = .95, col = "red")
  points(seq(1/6, 5/6, by = 1/6), data[ind, 1:5], col = c(3:6, 8), pch = plot.char, cex = 1.2)
  ind <- ind + 1
  for(i in 2:sim.per.row){
    
    plot.char <- rep(19, 5)
    plot.char[which.min(abs(data[ind, 1:5] - .95))] <- 1  
    
    plot(-1, -1, xlim = c(0, 1), ylim = c(.7, 1), xlab = "" ,ylab = "",
         xaxt = "n", yaxt = "n")
    abline(h = seq(.5, 1, by = .05), lty = 2, col = "gray", lwd = .5)
    abline(h = .95, col = "red")
    points(seq(1/6, 5/6, by = 1/6), data[ind, 1:5], col = c(3:6, 8), pch = plot.char, cex = 1.2)
    mtext(text = paste((data[ind, 21:25]), collapse = "."), line = 0, cex = .5)
    ind <- ind + 1
  }
}

legend.names <- paste(c("BLPR", "HDI", "SSLASSO", "SILM", "RR"), " (", win.rateA , ", ", win.rateN, ") ", sep = "")

reset()
par(xpd = TRUE)
legend("bottom", col = c(3:6, 8), pch = rep(19, 5),
       legend = legend.names, xpd =T, ncol = 5, cex = .8)

dev.off()

