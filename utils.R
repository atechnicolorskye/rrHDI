rm(list=ls())
pacman::p_load(tidyverse, xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

# Compile results into single csv file
setwd('../rr_perm_med_fixed')

dirnames <- list.dirs()
dirnames <- dirnames[-1]
for (d in dirnames){
  filenames <- list.files(path = d, full.names = TRUE)
  df <- filenames %>% lapply(read_csv) %>% bind_rows
  write.csv(df[-1], paste(substr(d, 3, nchar(d)), ".csv", sep = ""))
}

# experiment <- '100_300'
experiment <- 'rr_sign_50_100'
# filenames <- list.files(pattern = "_sign_", full.names = TRUE)
filenames <- list.files(full.names = TRUE)
df_means <- data.frame()
df_quantiles <- data.frame()
for (f in filenames){
  file <- read_csv(f)
  # browser()
  m_row <- cbind(t(colMeans(file[, 2:49])), file[1, 50:59])
  # med
  # m_row <- cbind(experiment, file[1, c(50, 51, 53)], t(colMeans(file[, c(13, 25, 37, 49)])))
  # small
  # m_row <- cbind(experiment, file[1, c(50, 51, 53)], t(colMeans(file[, c(9, 21, 33, 45)])))
  # colnames(m_row) <- c('Type', colnames(m_row)[2:4], substr(colnames(m_row)[5:8], 1, 5))
  # m_row <- cbind(experiment, file[1, c(6, 7, 9)], t(colMeans(file[, 2:5])))
  # colnames(m_row) <- c('Type', colnames(m_row)[2:4], substr(colnames(m_row)[5:8], 1, 5))
  # rows <- data.frame()
  # for (i in 1:3){
  #   m_row <- cbind(paste(substr(colnames(file[1, (1+i)]), 7, nchar(colnames(file[1, (1+i)]))), 'perm', experiment, sep = "_"),
  #                  file[1, c(14, 15, 17)], t(colMeans(file[, c(1+i, 4+i, 7+i, 10+i)])))
  #   colnames(m_row) <- c('Type', colnames(m_row)[2:4], substr(colnames(m_row)[5:8], 1, 5))
  # rows <- rbind(rows, m_row)
  # }
  # m_row <- rows
  df_means <- rbind(df_means, m_row)
  q_rows <- cbind(
            file[1, c(50, 51, 53)],
            paste(experiment, substr(names(file)[c(14:25, 38:49)], 1, nchar(names(file)[c(14:25, 38:49)])), sep = "_"),
            t(sapply(c(14:25, 38:49),
            # paste(experiment, substr(names(file)[c(21, 41)], 1, nchar(names(file)[c(21, 41)])), sep = "_"),
            # t(sapply(c(21, 41),
            # file[1, c(6, 7, 9)],
            # paste(experiment, substr(names(file)[c(3, 5)], 1, 5), sep = "_"),
            # t(sapply(c(3, 5),
            # file[1, c(14, 15, 17)],
            # paste(substr(names(file)[c(5:7, 11:13)], 7, nchar(names(file)[c(5:7, 11:13)])),
            #       'perm',
            #       experiment,
            #       substr(names(file)[c(5:7, 11:13)], 1, 5),
            #       sep = "_"),
            # t(sapply(c(5:7, 11:13),
            function(i){quantile(t(file[, i]), c(0.01, 0.25, 0.5, 0.75, 0.99))}))
           )
  colnames(q_rows) <- c(names(file)[c(50, 51, 53)], "Type", 0.01, 0.25, 0.5, 0.75, 0.99)
  # colnames(q_rows) <- c(names(file)[c(6, 7, 9)], "Type", 0.01, 0.25, 0.5, 0.75, 0.99)
  # colnames(q_rows) <- c(names(file)[c(14, 15, 17)], "Type", 0.01, 0.25, 0.5, 0.75, 0.99)
  df_quantiles <- rbind(df_quantiles, q_rows)
}
write.csv(df_means, "../rr_perm_med_all.csv")
# write.csv(df_quantiles, "../silm_perm_small_st_scale_ci_quantiles.csv")

filenames <- c("rr_perm_small_fixed.csv", "rr_sign_small_fixed.csv", "hdi_perm_small.csv", "hdi_sign_small.csv", "m_rr_hdi_perm_small.csv", "m_rr_hdi_sign_small.csv",
               "rr_perm_med_fixed.csv", "rr_sign_med_fixed.csv", "hdi_perm_med.csv", "hdi_sign_med.csv", "m_rr_hdi_perm_med.csv", "m_rr_hdi_sign_med.csv")
# filenames <- c("rr_perm_small_fixed_ci_quantiles.csv", "rr_sign_small_fixed_ci_quantiles.csv", "hdi_perm_small_ci_quantiles.csv", "hdi_sign_small_ci_quantiles.csv", 
#                "m_rr_hdi_perm_small_ci_quantiles.csv", "m_rr_hdi_sign_small_ci_quantiles.csv", 
#                "rr_perm_med_fixed_ci_quantiles.csv", "rr_sign_med_fixed_ci_quantiles.csv", "hdi_perm_med_ci_quantiles.csv", "hdi_sign_med_ci_quantiles.csv", 
#                "m_rr_hdi_perm_med_ci_quantiles.csv", "m_rr_hdi_sign_med_ci_quantiles.csv")
# filenames <- c("rr_perm_small_hard.csv", "hdi_perm_small_hard.csv", "m_rr_hdi_perm_small_hard.csv",
#                "rr_perm_med_hard.csv", "hdi_perm_med_hard.csv", "m_rr_hdi_perm_med_hard.csv")
# filenames <- c("rr_perm_small_hard_ci_quantiles.csv", "hdi_perm_small_hard_ci_quantiles.csv", "m_rr_hdi_perm_small_hard_ci_quantiles.csv",
#                "rr_perm_med_hard_ci_quantiles.csv", "hdi_perm_med_hard_ci_quantiles.csv", "m_rr_hdi_perm_med_hard_ci_quantiles.csv")

combined_means <- data.frame()
for (f in filenames){
  file <- read_csv(f)
  combined_means <- rbind(combined_means, file)
}
write.csv(combined_means[-1], "combined_means.csv", row.names = FALSE)
# write.csv(combined_means[-1], "combined_quantiles_hard.csv", row.names = FALSE)


# df <- filenames %>% lapply(read_csv) %>% bind_rows
# df <- df[, -1]
# write.csv(df, '../rr_perm_small_fixed.csv')
# 
# filenames <- list.files(pattern = "n_100_", full.names = TRUE)
# df <- filenames %>% lapply(read_csv) %>%  bind_rows
# df <- df[, -1]
# write.csv(df, 'large.csv')

## Process results
# *_admissible are results filtered for admissible group actions
# *_admissible_drop are results without unscaled RR
# data <- read.csv("gaussian_jm.csv", stringsAsFactors = F)
excel_sheets("partial_perm.xlsx")
data <- read_excel("partial_perm.xlsx", sheet = "small_perm_active")
colnames(data) <- c("s", "X", "e", "A C BLPR", "A C HDI", "A C DL", "A C SILM", "A C RR", "A L BLPR", "A L HDI", "A L DL", "A L SILM", "A L RR", "I C BLPR", "I C HDI", "I C DL", "I C SILM", "I C RR", "I L BLPR", "I L HDI", "I L DL", "I L SILM", "I L RR")
print(xtable(data, align=rep('c', 24), digits=c(0, 0, 0, 0, rep(3, 20))), include.rownames=F)

# data <- read.csv("out/b_1000_n_solve_500/large_admissible_drop.csv", stringsAsFactors = F)
# data <- read.csv("out/b_1000_n_solve_500/small_admissible_drop.csv", stringsAsFactors = F)
data[, 1:20] <- as.numeric(as.matrix(data[, 1:20]))
win.rateA <- tabulate(apply(data[, 1:5], MAR = 1, function(x){which.min(abs(x - .95))}), nbins=5)
win.rateN <- tabulate(apply(data[, 11:15], MAR = 1, function(x){which.min(abs(x - .95))}), nbins=5)

# Min, Q1, Q2, Q3, Max
tab1 <- round(cbind(apply(data[, 1:5], MAR = 2, min, na.rm =T),
                    apply(data[, 1:5], MAR = 2, function(x){quantile(x, .25)}),
                    apply(data[, 1:5], MAR = 2, median, na.rm =T),
                    apply(data[, 1:5], MAR = 2, function(x){quantile(x, .75)}),
                    apply(data[, 1:5], MAR = 2, max, na.rm =T),
                    apply(data[, 11:15], MAR = 2, min, na.rm =T),
                    apply(data[, 11:15], MAR = 2, function(x){quantile(x, .25)}),
                    apply(data[, 11:15], MAR = 2, median, na.rm =T),
                    apply(data[, 11:15], MAR = 2, function(x){quantile(x, .75)}),
                    apply(data[, 11:15], MAR = 2, max, na.rm =T)), 3)

colnames(tab1) <- c("min CR-A", "Q1 CR-A", "Q2 CR-A", "Q3 CR-A", "Max CR-A", "min CR-N", "Q1 CR-N", "Q2 CR-N", "Q3 CR-N", "Max CR-N")
rownames(tab1) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")

print(xtable(tab1, digits=3, include.rownames=F))

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
rownames(tab2) <- c("BLPR", "HDI", "SSLASSO", "SILM", "RR")

# Rankings
rank_dist_a <- apply(data[, 1:5], MAR = 1, function(x){round(abs(x - .95), 5)})
tab3 <- data.frame(data$s, data$x, data$e, data$g, t(apply(rank_dist_a, MAR = 2, rank)))
colnames(tab3) <- c("s", "Dist of X", "Dist of eps", "Invariance", "BLPR", "HDI", "SSLASSO", "SILM", "RR")
means <- round(colMeans(t(apply(rank_dist_a, MAR = 2, rank))), 2)
# tab3 <- rbind(tab3, c(rep("", 4), means))
# rownames(tab3)[rownames(tab3) == "means"] <- "Average Rank for Active Parameters"

xtable(matrix(means, nrow = 1))
print(xtable(tab3), include.rownames = F)

rank_dist_n <- apply(data[, 11:15], MAR = 1, function(x){round(abs(x - .95), 5)})
tab4 <- data.frame(data$s, data$x, data$e, data$g, t(apply(rank_dist_n, MAR = 2, rank)))
colnames(tab4) <- c("s", "Dist of X", "Dist of eps", "Invariance", "BLPR", "HDI", "SSLASSO", "SILM", "RR")
means <- round(colMeans(t(apply(rank_dist_n, MAR = 2, rank))), 2)

# tab4 <- rbind(tab4, means)
# rownames(tab4)[rownames(tab4) == "means"] <- "Average Rank for Inactive Parameters"

xtable(matrix(means, nrow = 1))
print(xtable(tab4), include.rownames = F)

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

