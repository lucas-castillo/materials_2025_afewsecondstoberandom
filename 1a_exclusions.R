rm(list=ls())
library(tidyverse)
library(glue)
source("src/RQA.R")
source("src/exclusion_utils.R")
# Load data
data <- read_csv(glue("data/speed/rg_table.csv"))

# Add determinism as exclusion criterion (DET > .75) in either DET of seq or DET of diff(seq)
data <- data %>% 
  group_by(uuid, condition) %>% 
  # set determinism
  mutate(DET = max(RQA(value)$DET, RQA(diff(value))$DET)) %>% 
  mutate(tooPredictable = DET > .75)

# Remove problematic values: lifespans above 140 or 4SD > mean
data <- data %>% 
  group_by(uuid) %>% 
  mutate(ub = mean(value, na.rm = T) + 4 * sd(value, na.rm=T)) %>% 
  mutate(ub = ifelse(ub > 140, 140, ub)) %>% # never allow upper bound above 140
  ungroup() %>% 
  mutate(value2 = ifelse(value > ub, NA, value)) %>% 
  mutate(value = value2) %>% 
  select(-c(ub, value2))  
# Add out of range as exclusion criterion in speed
hard_lb <- 20
hard_ub <- 200
soft_lb <- -1
soft_ub <- 130
threshold <- .05

data <- data %>% 
  group_by(uuid, condition) %>% 
  mutate(outOfRange = out_of_range(
    value, threshold=threshold,
    hard_lb=hard_lb, hard_ub=hard_ub, 
    soft_lb=soft_lb, soft_ub=soft_ub
  )
)
# Add production speed as exclusion criterion for speed
# (see slow_pace_harsh for original, preregistered criterion)
data <- data %>% 
  mutate(expect = ifelse(condition == "Fast", 80, 40)) %>% 
  group_by(uuid, condition) %>% 
  mutate(good_window = slow_pace_p(value, expect*2, .6)) %>% 
  mutate(tooSparse = sum(good_window) == 0) %>% 
  select(-expect)
data %>% filter(!outOfRange & !tooPredictable & !tooSparse) %>% 
  select(-c(outOfRange, tooSparse, tooPredictable, DET, good_window)) %>% 
  ungroup %>% 
  save(file = glue("data/rdata/speed.RData"))
