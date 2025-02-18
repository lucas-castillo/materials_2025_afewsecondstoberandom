rm(list = ls())
library(tidyverse)
library(glue)
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
all <- tibble()
bpms <- tibble()
for (exp in exps){
  data <- get(load(glue("data/rdata/{exp}.RData")))
  data <- data %>% 
    distinct(uuid, condition, block, bpm) %>% mutate(exp=exp)
  data <- data %>% mutate(across(uuid:block, as.character))
  bpms <- bpms %>% rbind(data)
}

for (exp in exps){
  print(exp)
  fpath <- glue("output/ratio-measures/{exp}/processed/eps")
  ind_files <- list.files(fpath)
  ind_files <- ind_files[ind_files != "normalizing_values.RData" & ind_files != "eps"]
  last_increment <- -1
  for (f in ind_files){
    increment <- which(f == ind_files) / length(ind_files) * 100
    if (increment > (last_increment + 5)){
      cat(increment)
      cat("||")
      last_increment <- increment
    }
    
    search_bpm <- function(u,c,b){
      if (is.na(b) | b == "NA"){
        return(bpms$bpm[bpms$uuid == u & bpms$condition == c])
      } else{
        return(bpms$bpm[
          bpms$uuid == u & 
            bpms$condition == c &
            bpms$block == b
          
        ])
        
      }
    }
    
    load(glue("{fpath}/{f}"))
    all <- all %>%
      rbind(
        temp_person.p %>%
          ungroup %>% 
          mutate(exp=exp) %>%
          nest(.by=c(uuid, condition, block, exp)) %>%
          mutate(bpm = search_bpm(uuid, condition, block)) %>%
          unnest(data)
      )
  }
  cat("\n---\n")
}

save(all, file = "output/ratio-measures/people_p.RData")

