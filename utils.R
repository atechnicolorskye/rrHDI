rm(list=ls())
pacman::p_load(tidyverse, xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


##  Compile results into single csv file
# set working directory
setwd('~/out')
# comb through directories, combines files and cleans up
dirnames_method <- list.dirs(path = '.', recursive = FALSE)
for (dm in dirnames_method){
  dirnames <- list.dirs(path = dm, full.names = TRUE, recursive = FALSE)
  for (d in dirnames){
    filenames <- list.files(path = d, full.names = TRUE)
    df <- filenames %>% lapply(read_csv) %>% bind_rows
    write.csv(df[-1], paste(substr(d, 3, nchar(d)), ".csv", sep = ""))
    unlink(d, recursive = TRUE, force = TRUE)
  }  
}

## Output summary statistics used to generate plots
dirnames_method <- list.dirs(path = '.', recursive = FALSE)
for (dm in dirnames_method){
  filenames <- list.files(path = dm, full.names = TRUE)
  df_means <- data.frame()
  df_quantiles <- data.frame()
  for (f in filenames){
    file <- read_csv(f)
    m_row <- cbind(t(colMeans(file[, c(2:13)])), substr(names(file)[2], 9, nchar(names(file)[2])), file[1, c(14:15, 17, 20:22)])
    colnames(m_row)[1:13] <- c(substr(colnames(m_row)[1:12], 1, 7), "Method")
    df_means <- rbind(df_means, m_row)
    q_rows <- cbind(
      file[1, c(14, 15, 17, 22)],
      substr(names(file)[c(5:7, 11:13)], 1, 7),
      substr(names(file)[c(5:7, 11:13)], 9, nchar(names(file)[c(5:7, 11:13)])),
      t(sapply(c(5:7, 11:13),
               function(i){quantile(t(file[, i]), c(0.01, 0.25, 0.5, 0.75, 0.99))}))
    )
    colnames(q_rows) <- c(names(file)[c(14, 15, 17, 21)], "Type" ,"Method", 0.01, 0.25, 0.5, 0.75, 0.99)
    df_quantiles <- rbind(df_quantiles, q_rows)
  }
  write.csv(df_means, paste(substr(dm, 3, nchar(dm)), ".csv", sep=""))
  write.csv(df_quantiles, paste(substr(dm, 3, nchar(dm)), "_ci_quantiles.csv", sep=""))
}  

## Output tables 
# get summary csvs
files <- list.files(path='.', pattern='.csv', recursive = FALSE)
files_ci <- list.files(path='.', pattern='_ci_quantiles.csv', recursive = FALSE)
files_cov <- setdiff(files, files_ci)

# set s
s <- 4
# set quantile of interest from [0.25, 0.5, 0.75, 0.99]
ci_q <- '0.75'

# coverage
combined <- data.frame()
for (f in 1:length(files_cov)){
  file <- read_csv(files_cov[f])
  file$X1 <- substr(files_cov[f], 1, nchar(files_cov[f]) - 4)
  combined <- rbind(combined, file)
}
combined <- combined[which(combined$s == s), ]
combined <- combined[with(combined, order(combined$x, combined$e)), ]
combined <- combined[c('X1', 'a_cov_i', 'a_cov_a', 'a_cov_s', 'n_cov_i', 'n_cov_a', 'n_cov_s', 'x', 'e')]
print(xtable(combined %>% mutate(across(where(is.numeric), round, 2))), include.rownames=F)

# confidence interval quantiles
combined <- data.frame()
for (f in 1:length(files_ci)){
  file <- read_csv(files_ci[f])
  file$X1 <- substr(files_ci[f], 1, nchar(files_ci[f]) - 17)
  combined <- rbind(combined, file)
}
combined <- combined[which(combined$s == s), ]
combined <- combined[with(combined, order(combined$x, combined$e)), ]
combined <- combined[c('X1', 'e', 'Type', ci_q)]
combined_n <- data.frame()
for (i in unique(combined$Type)){
  if (nrow(combined_n) == 0 && ncol(combined_n) == 0){
    combined_n <- combined[which(combined$Type == i), ][c('X1', ci_q)]
  }
  else{
    combined_n <- cbind(combined_n, combined[which(combined$Type == i), ][c(ci_q)])
  }
}
colnames(combined_n) <- c(colnames(combined_n)[1], unique(combined$Type))
print(xtable(combined_n %>% mutate(across(where(is.numeric), round, 2))), include.rownames=F)
