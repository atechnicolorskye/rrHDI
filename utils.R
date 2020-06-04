rm(list=ls())
library(tidyverse)
# setwd('out')

filenames <- list.files(pattern = "_n_50.+perm", full.names = TRUE)
df <- filenames %>% lapply(read_csv) %>%  bind_rows
df <- df[, -1]
write.csv(df, 'small_perm_no_scale.csv')

filenames <- list.files(pattern = "sign.+_n_100_", full.names = TRUE)
df <- filenames %>% lapply(read_csv) %>%  bind_rows
df <- df[, -1]
write.csv(df, 'large.csv')

# filenames <- list.files(pattern = "_n_500_.+_i_-", full.names = TRUE)
# df <- filenames %>% lapply(read_csv) %>%  bind_rows
# df <- df[, -1]
# write.csv(df, 'small.csv')

