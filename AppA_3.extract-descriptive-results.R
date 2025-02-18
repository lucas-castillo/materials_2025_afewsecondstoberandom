rm(list=ls())
library(tidyverse)
library(glue)
library(lmerTest)
library(rstanarm)
library(logspline)
library(bayestestR)
set.seed(2023)

data <- get(load(glue("data/rdata/fs_fulldescriptives.RData"))) 
data_sv <- get(load(glue("data/rdata/fs_sv.RData"))) 
measure_names <- data %>% 
  select(contains("boot_")) %>% 
  colnames %>% 
  str_remove("boot_")  
logit_measures <- c("R", "A", "TPF")
sv_measures <- data_sv %>% 
  select(-c(uuid:blockorder), -contains("boot")) %>% 
  colnames

print_results <- function(m, expectations=F){
  f_model <- glue("mod{m}")
  bf <- glue("mod{m}.bf")
  is_logit <- m %in% logit_measures
  if (expectations | !is_logit){
    cnames <- c("beta", "SE", "df", "t", "p")
  } else {
    cnames <- c("beta", "SE", "Z", "p")
  }

  f_table <- get(f_model) %>% 
    summary %>% 
    coefficients %>% 
    round(5) %>% 
    magrittr::set_colnames(cnames)
  b_table <- get(bf)$log_BF %>% 
    effectsize::interpret_bf(log=T, include_value = T) %>% 
    data.frame(bf = .)
    
  cbind(
    f_table,
    b_table
  ) %>% 
  knitr::kable()
}

exps <- c("fs")
for (exp in exps){
  print(exp)
  path <- glue("output/descriptive-analyses/{exp}")
  analyses <- str_remove(list.files(path)[str_detect(list.files(path), ".RData")], ".RData")
  for (analysis in analyses){
    print(analysis)
    lf <- load(glue("{path}/{analysis}.RData"))
    write("", glue("{path}/{analysis}.txt"), append = F)
    for (m in measure_names){
      
      print(m)
      write(glue("## {m} --------------"), glue("{path}/{analysis}.txt"), append = T)
      write(print_results(m, expectations = analysis=="expectations-condition"), glue("{path}/{analysis}.txt"), append = T)
      write("\n", glue("{path}/{analysis}.txt"), append = T)  
    }
    if (analysis == first(analyses)){
      data <- modR.b$data
      write("**Block**\n", glue("{path}/contrasts.txt"), append = F)
      write(knitr::kable(contrasts(data$block)), glue("{path}/contrasts.txt"), append = T)
      write("\n**Block Order**\n", glue("{path}/contrasts.txt"), append = T)
      write(knitr::kable(contrasts(data$blockorder)), glue("{path}/contrasts.txt"), append = T)
      write("\n**Condition**\n", glue("{path}/contrasts.txt"), append = T)
      write(knitr::kable(contrasts(data$condition)), glue("{path}/contrasts.txt"), append = T)
      write("\n", glue("{path}/contrasts.txt"), append = T)
    }
    rm(list = c(lf, "lf"))
  }
}


