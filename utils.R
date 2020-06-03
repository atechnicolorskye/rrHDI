setwd('out')

filenames <- list.files(pattern = "_n_50_", full.names = TRUE)
df <- filenames %>% lapply(read_csv) %>%  bind_rows
df <- df[, -1]
write.csv(df, 'small.csv')

filenames <- list.files(pattern = "_n_100_", full.names = TRUE)
df <- filenames %>% lapply(read_csv) %>%  bind_rows
df <- df[, -1]
write.csv(df, 'large.csv')