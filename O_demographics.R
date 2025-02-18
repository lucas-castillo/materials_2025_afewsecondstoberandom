rm(list=ls())
library(tidyverse)
d <- read_csv("data/speed/survey_table.csv")
df <- get(load("data/rdata/speed.RData"))

# Age
d %>% 
  filter(uuid %in% unique(df$uuid)) %>% 
  mutate(age = as.numeric(ifelse(age == "Prefer not to say", NA, age))) %>% 
  summarise(M = mean(age, na.rm = T), SD = sd(age, na.rm = T))
# Gender
d %>% 
  filter(uuid %in% unique(df$uuid)) %>% 
  group_by(gender) %>% 
  tally
## % participants with both sequences
df %>% 
  ungroup %>% 
  distinct(uuid, condition) %>% 
  group_by(uuid) %>% 
  tally() %>% 
  group_by(n) %>% 
  tally %>% 
  mutate(snn = sum(nn)) %>% 
  mutate(p = nn / snn)  
