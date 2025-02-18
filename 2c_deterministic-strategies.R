rm(list=ls())
library(tidyverse)
library(glue)
set.seed(2023)
source("src/measures_after_thinning.R")
# get Ns from data --------------------------------------------------------
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
  "groot22"
)
Ns <- tibble()
for (exp in exps){
  data <- get(load(glue("data/rdata/{exp}.RData")))
  N <- data %>% 
    group_by(uuid, condition, block) %>% 
    tally %>% 
    ungroup %>% 
    summarise(n = mean(n)) %>% 
    pull(n)
  
  if (str_detect(exp, "castillo24_e1|speed")) {
    domain <- "naturalistic"
  } else if (str_detect(exp, "cooper|towse|castillo24_e2")) {
    domain <- "range"
  } else{
    domain <- "binary"
  }
  Ns <- rbind(
    Ns, 
    tibble(domain, N)
  )
}
Ns <- Ns %>% 
  group_by(domain) %>% 
  summarise(N = round(mean(N)))

# set deterministic strategies --------------------------------------------
domains <- Ns$domain
strategies <- c("sweep", "sweepback", "rep_rand", "edge_center")
generate_det <- function(
    domain = "range",
    strategy = "",
    N = 100,
    NN=10
){
  if (strategy == "sweep"){
    ## ABCDEFGABCDEFG
    vv <- unique(sort(generate_random(domain, N = NN, interval = 1)))
    return(rep(vv, ceiling(N / length(vv)))[1:N])
  } else if (strategy == "sweepback"){
    ## ABCDEFGFEDCBA
    vv <- unique(sort(generate_random(domain, N = NN, interval = 1)))
    return(rep(c(vv, sort(vv[2:(length(vv) - 1)], T)), ceiling(N / length(vv)))[1:N])
  } else if (strategy == "rep_rand"){
    ## repeat a random string: DAEFCGB,DAEFCGB,DAEFCGB
    vv <- unique(generate_random(domain, N = NN, interval = 1))
    return(rep(vv, ceiling(N / length(vv)))[1:N])
  } else if (strategy == "edge_center"){
    ## AGBFCED
    vv <- unique(sort(generate_random(domain, N = NN, interval = 1)))
    new_vv <- c()
    for (i in 1:length(vv)){
      new_vv[i] <- if (i %% 2 == 1) vv[floor(i/2) + 1] else vv[length(vv) - floor(i/2) + 1]
    }
    return(rep(new_vv, ceiling(N / length(vv)))[1:N])
  }
}

for (domain in domains){
  data <- tibble()
  for (strategy in strategies){
    for (NN in c(10, 20)){ # strategy horizon
      N <- Ns$N[Ns$domain == domain]
      for (sim in 1:1000){
        data <- rbind(
          data, 
          tibble(
            index = 1:N,
            uuid=glue("det.{domain}.{strategy}.{NN}.{sim}"),
            condition=glue("det_strategies.{domain}.{strategy}.{NN}"),
            block=1,
            value=generate_det(domain, strategy, N, NN),
            bpm=100
          )
        )
      }
    }
  }
  save(data, file = glue("data/rdata/det_{domain}.RData"))
}


