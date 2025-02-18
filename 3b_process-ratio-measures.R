rm(list = ls())
library(tidyverse)
library(glue)
library(lmerTest)
library(rstanarm)
options(mc.cores = parallel::detectCores())
theme_set(theme_minimal(15))
source("src/measures_after_thinning.R")
measure_names <- c(
  "R", 
  "A",
  "TPF", 
  "D", 
  "RNG", 
  "SWLZ_ER", 
  "gzip", 
  "lz76"
)
rerun <- T
exps <- c(
  "cooper12",
  "towse16_e1",
  "towse16_e2",
  "towse16_e3",
  "castillo24_e1",
  "castillo24_e2",
  "speed2",
  "groot22",
  "guseva23_er",
  "guseva23_irr",
  "guseva23_mt",
  "half_speed",
  "det_binary",
  "det_range",
  "det_naturalistic"
)
safe_sep <- function(x, into, sep){
  split <- x %>% 
    stringi::stri_reverse() %>% # reverse the string
    str_split_1(sep) %>% # split it
    sapply(stringi::stri_reverse, USE.NAMES = F) # reverse all back
  split <- rev(split)
  tibble(v=split) %>% 
    mutate(names = into) %>% 
    pivot_wider(names_from = names, values_from = v)
  
  # t(bind_cols(sapply(str_split_1(stringi::stri_reverse(x), "__"), stringi::stri_reverse, USE.NAMES = F)))
}

for (exp in exps){
  print(exp)
  dir.create(glue("output/ratio-measures/{exp}/processed/eps/"), showWarnings = F, recursive = T)
  isFirst <- T
  sequence_files <- list.files(glue("output/ratio-measures/{exp}"))
  sequence_files <- sequence_files[str_detect(sequence_files, ".csv")]
  print(length(sequence_files))
  # Get M+SD of observed values, used to normalize everything later ----------
  if (isFirst){ # after the first file this can be skipped as already in environment
    normalizing_values <- tibble()
    for (f in sequence_files){
      temp <- read_csv(glue("output/ratio-measures/{exp}/{f}"), show_col_types = F) %>%
        filter(source == "participant" & type == "full") %>%
        mutate(sequence_id = str_remove(f, ".csv")) %>% 
        mutate(safe_sep(sequence_id, into=c("uuid", "condition", "block"), sep ="__"))
      normalizing_values <- rbind(normalizing_values, temp)
    }
    normalizing_values <- normalizing_values %>%
      group_by(condition) %>%
      reframe(across(all_of(measure_names), \(x){c("M"=mean(x, na.rm=T), "SD"=sd(x, na.rm=T))})) %>%
      mutate(summary=rep(c("M", "SD"), length(unique(normalizing_values$condition)))) %>%
      pivot_longer(all_of(measure_names), names_to = "measure")
  }
  # Process simulation, one uuid-condition-block at a time -------------------
  for (f in sequence_files){
    print(glue(
      "{which(f == sequence_files)} / {length(sequence_files)}"
    ))
    
    info <- f %>% 
      str_remove(".csv") %>% 
      safe_sep(into = c("uuid", "condition", "block"), sep = "__")
    uuid <- info$uuid
    condition <- info$condition
    block <- info$block
    
    M <- read_csv(glue("output/ratio-measures/{exp}/{f}"), show_col_types = F)
    result <- process_ratio_measures(
      M, measure_names, uuid, condition, block, exp=exp
    )
    
    temp_iid <- result[[1]]
    temp_iid.distance <- result[[2]]
    temp_iid_summaries <- result[[3]]
    temp_person <- result[[4]]
    temp_person.p <- result[[5]]
    
    save(
      temp_iid, temp_iid.distance, temp_iid_summaries, temp_person, temp_person.p,
      file=glue("output/ratio-measures/{exp}/processed/eps/{str_remove(f, '.csv')}.RData")
    )
    if (isFirst){
      isFirst <- !isFirst
      save(
        normalizing_values,
        file=glue("output/ratio-measures/{exp}/processed/eps/normalizing_values.RData")
      )
    }
  }
}
