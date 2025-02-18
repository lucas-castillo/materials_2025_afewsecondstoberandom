library(tidyverse)
library(glue)
# Setup -------------------------------------------------------------------
source("src/get_1D_summaries.R")
set.seed(2023)
n_bootstrap <- 1e4
# Load Data
data <- get(load(glue("data/rdata/speed.RData")))

bootstrap_exp <- function(seq){
  # Generates tibble of expectations for each measure
  M <- tibble()
  a <- min(seq[!is.na(seq)]):max(seq[!is.na(seq)])
  for (i in 1:n_bootstrap){
    temp_seq <- sample(seq)
    M <- rbind(M, all_measures(temp_seq, a))
  }
  M %>% 
    summarise(across(everything(), \(x){mean(x, na.rm = T)})) %>% 
    pivot_longer(everything()) %>% 
    mutate(name = str_c("boot_", name)) %>% 
    pivot_wider(names_from = name, values_from = value)
}

get_long_values <- function(seq){
  # makes df with values that can be described on a per-item level (e.g. repetitions)
  a <- min(seq[!is.na(seq)]):max(seq[!is.na(seq)])
  all_measures_full(seq, a = a)[[1]]
}
get_short_values <- function(seq){
  # makes df with values that can only be described on a per-sequence level (e.g. compression measures)
  a <- min(seq[!is.na(seq)]):max(seq[!is.na(seq)])
  all_measures_full(seq, a = a)[[2]]
}

# Change variables to factors ---------------
data <- data %>% 
  mutate(
    uuid = factor(uuid),
    condition = factor(condition),
    block = factor(block),
    blockorder = factor(blockorder)
  )

# Add Bootstrap Expectations -------
data <- data %>% 
  group_by(uuid, condition) %>% 
  mutate(bind_cols(bootstrap_exp(value)))

# Add long values ---------
data <- data %>% 
  group_by(uuid, condition) %>% 
  mutate(bind_cols(get_long_values(value))) %>% 
  ungroup

# Make a df for short values --------------
# I'm sure there's a clever way to get the boot_* colnames programatically, but they're hardcoded here
data_sv <- data %>% 
  group_by(across(c(uuid, condition, block, blockorder, boot_R:boot_lz76))) %>% 
  do(
    get_short_values(.$value)
  ) %>% 
  ungroup

# Do the same for a version of the data where the fast sequences get truncated to have same
# number of items as slow sequences (as some measures are length-dependent)
two_seqs_uuids <- data %>% 
  distinct(uuid, condition) %>% 
  filter(condition %in% c("Fast", "Slow")) %>% 
  group_by(uuid) %>% 
  tally %>% 
  filter(n == 2) %>% 
  pull(uuid)

shorter_ns <- data %>% 
  filter(condition %in% c("Fast", "Slow")) %>% 
  filter(uuid %in% two_seqs_uuids) %>% 
  group_by(uuid, condition) %>% 
  tally %>% 
  group_by(uuid) %>% 
  mutate(shorter_n = min(n)) %>% 
  distinct(uuid, shorter_n)

data.t <- data %>% 
  # truncate sequences --
  filter(condition %in% c("Fast", "Slow")) %>% 
  filter(uuid %in% two_seqs_uuids) %>% 
  nest(.by=uuid) %>% 
  rowwise() %>% 
  mutate(target_n = shorter_ns$shorter_n[shorter_ns$uuid == uuid]) %>% 
  unnest(data) %>% 
  filter(index <= target_n)

data_sv.t <- data.t %>% 
  group_by(across(c(uuid, condition, block, blockorder, boot_R:boot_lz76))) %>% 
  do(
    get_short_values(.$value)
  ) %>% 
  ungroup



data <- data %>% ungroup()
data_sv <- data_sv %>% ungroup()
data.t <- data.t %>% ungroup() %>% 
    mutate(condition = factor(condition))
data_sv.t <- data_sv.t %>% ungroup() %>% 
  mutate(condition = factor(condition))

# Fast vs Slow
data.fs <- data %>% 
  filter(condition %in% c("Fast", "Slow")) %>% 
  mutate(condition = factor(condition))
data_sv.fs <- data_sv %>% 
  filter(condition %in% c("Fast", "Slow")) %>% 
  mutate(condition = factor(condition))

comparison <- "fs"
save(list=glue("data.{comparison}"), file = glue("data/rdata/{comparison}_fulldescriptives.RData"))
save(list=glue("data_sv.{comparison}"), file = glue("data/rdata/{comparison}_sv.RData"))
if (comparison == "fs"){
  save(list=glue("data.t"), file = glue("data/rdata/fs_fulldescriptives_trunc.RData"))
  save(list=glue("data_sv.t"), file = glue("data/rdata/fs_sv_trunc.RData"))
}  
