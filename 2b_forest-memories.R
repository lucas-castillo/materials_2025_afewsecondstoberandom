rm(list=ls())
library(tidyverse)
library(glue)
library(abcrf)
library(ranger)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())
select <- dplyr::select

source("src/rforest_memories.R")

## produce leaf memories from forests (for posteriors) --------------
load(file = "output/ratio-measures/recovery/abcrf_models.RData")
mem_b <- get_forest_memory(
  model_b,
  training = sims %>% filter(condition == "binary")
)
saveRDS(mem_b, "output/ratio-measures/recovery/mem_b.rds")

mem_r <- get_forest_memory(
  model_r,
  training = sims %>% filter(condition == "range")
)
saveRDS(mem_r, "output/ratio-measures/recovery/mem_r.rds")

mem_n <- get_forest_memory(
  model_n,
  training = sims %>% filter(condition == "naturalistic")
)
saveRDS(mem_n, "output/ratio-measures/recovery/mem_n.rds")