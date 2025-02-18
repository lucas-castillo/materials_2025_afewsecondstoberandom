rm(list=ls())
library(randMeasures)
library(tidyverse)
library(glue)
library(posterior)
library(foreach)
library(doParallel)
source("src/measures_after_thinning.R")
registerDoParallel(parallel::detectCores())
set.seed(2023)
possible_ratios <- 1:7

exps <- c(
  "guseva23_er",
  "guseva23_irr",
  "guseva23_mt",
  "cooper12", 
  "towse16_e1",
  "towse16_e2",
  "towse16_e3",
  "castillo24_e1",
  "castillo24_e2",
  "speed2",
  "groot22",
  "half_speed",
  "det_binary",
  "det_range",
  "det_naturalistic"
)

# Compute measures after thinning -----------------------------------------
for (exp in exps){
  print(exp)
  dir.create(glue("output/ratio-measures/{exp}/"), showWarnings = F, recursive = T)
  data <- get(load(glue("data/rdata/{exp}.RData")))
  data <- data %>% 
    mutate(sequence_id=paste(uuid, condition, block, sep = "__"))
  data <- data %>% arrange(uuid, condition, block)
  is.circular <- exp == "cooper12"
  
  all <- foreach(
    seq_id=unique(data$sequence_id), .export = "exp",
    .packages = c("dplyr", "magrittr", "tidyr", "glue", "readr", "stringr", "randMeasures")
  ) %dopar% {
    temp <- data %>% filter(sequence_id == seq_id)
    ratio_measures <- compute_ratio_thinning(
      v = temp$value,
      ratios = possible_ratios,
      do_iid = T,
      N_IID = 2e3,
      is.circular = is.circular,
      drop.r = F
    )

    ratio_measures <- ratio_measures %>%
      mutate(bpm = temp$bpm[1])
    write_csv(ratio_measures, glue("output/ratio-measures/{exp}/{seq_id}.csv"))
  }
}
