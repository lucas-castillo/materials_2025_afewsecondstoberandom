rm(list=ls())
library(abcrf)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())
library(tidyverse)
library(glue)
select <- dplyr::select
source("src/measures_after_thinning.R")
dir.create(path = "output/ratio-measures/recovery/simulations/", recursive = T, showWarnings = F)
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

# Simulate sequences ------------------------------------------------------
set.seed(2024)
schemas <- 0:+5 # same in the negative so skip
dfuns <- c("dmean", str_c("d", schemas))

N_SIMS <- length(dfuns) * (length(dfuns) - 4) * length(8:9)
N_SIMS <- N_SIMS * 5
planned_simulations <- expand_grid(
  f_interval = c("NO", as.character(1:7), "H"),
  i = 1:2,
  dfun = dfuns
) %>% 
  # split f_interval=H into interval=8 and interval=9
  mutate(even = as.numeric(i %% 2 == 0)) %>% 
  mutate(interval = ifelse(f_interval == "NO", 0,
                    ifelse(f_interval == "H", 8 + even, f_interval))) %>% 
  mutate(interval = as.numeric(interval)) %>% 
  # remove impossible dfun 
  distinct(f_interval, interval, dfun) %>% 
  filter(!(interval == 0 & dfun %in% c("dmean", "d0", "d2", "d4"))) %>% 
  # calculate how many per cell so that all f_interval have the same
  group_by(f_interval) %>% 
  mutate(N = n()) %>% 
  mutate(n = round(N_SIMS / n()))
# Check: This must be true for a 1/N prior
planned_simulations %>% 
  group_by(f_interval) %>% 
  summarise(n = sum(n)) %>% 
  summarise(sd = sd(n)) %>% 
  pull(sd) %>% 
  magrittr::equals(0)

simulated_values <- tibble()
Ns <- Ns %>% mutate(N = round(N))
for (row_index in 1:nrow(planned_simulations)){
  dfun <- planned_simulations[[row_index, "dfun"]]
  interval <- planned_simulations[[row_index, "interval"]]
  n <- planned_simulations[[row_index, "n"]]
  for (domain in Ns$domain){
    if (dfun != "dmean"){
      num <- as.numeric(str_split_1(dfun, "d")[2])
      f <- \(x){x+num}
    } else {
      if (domain == "binary"){
        f <- \(x){1}  
      } else if (domain == "range"){
        f <- \(x){5}  
      } else if (domain == "naturalistic"){
        f <- \(x){170}  
      }
    }
    for (i in 1:n){
      simulated_values <- rbind(
        simulated_values,
        tibble(
          v = generate_random(
            interval=interval,
            N = Ns$N[Ns$domain == domain],
            domain = domain,
            det_f = f
          ),
          idx = 1:Ns$N[Ns$domain == domain],
          interval = interval,
          dfun=dfun,
          id = paste0(domain, "_", i, "_", interval, "_", dfun)
        )
      )
    }
  }
}

saveRDS(object=simulated_values, file="output/ratio-measures/recovery/simulated_values.rds")

# Compute measures --------------------------------------------------------
set.seed(2024)
simulated_values <- readRDS(file="output/ratio-measures/recovery/simulated_values.rds")
results <- foreach(
  id=unique(simulated_values$id),    
  .packages = c("dplyr", "magrittr", "tidyr", "glue", "readr", "stringr")
  ) %dopar% {
  i <- which(unique(simulated_values$id) == id)
  v <- simulated_values$v[simulated_values$id == id]
  M <- compute_ratio_thinning(
    v, ratios = 1:7,       
    do_iid = T,
    N_IID = 2e3,
    is.circular = F,
    drop.r = F
  )
  write_csv(M, glue("output/ratio-measures/recovery/simulations/{id}.csv"))
}

# Get empirical ps from measures ----------------------------------------------
set.seed(2024)
EPs <- foreach(
  f=list.files("output/ratio-measures/recovery/simulations/"), 
  .combine = "rbind",
  .packages = c("dplyr", "magrittr", "tidyr", "glue", "readr", "stringr", "randMeasures")
) %dopar% {
  sequence_id <- str_remove(f, ".csv")
  info <- str_split_1(sequence_id, pattern = "_")
  domain <- info[1]
  sim <- info[2]
  interval <- info[3]
  det_function <- info[4]
  process_ratio_measures(
    read_csv(glue("output/ratio-measures/recovery/simulations/{f}"), show_col_types = F) %>% 
      mutate(condition = domain), 
    measure_names, 
    uuid = sequence_id, condition = domain, block = 1, exp = domain
  )[[5]]
}
save(EPs, file = "output/ratio-measures/recovery/method_recovery_EPs.RData")

# Run ABCRF ----------------------------------------------------------
set.seed(2024)
load(file = "output/ratio-measures/recovery/method_recovery_EPs.RData")
EPs <- EPs %>% 
  ungroup %>% 
  nest(.by=uuid) %>% 
  separate(uuid, into = c("domain", "sim", "interval", "dfun"), sep="_", remove=F) %>% 
  unnest(data)


sims <- EPs %>% 
  group_by(across(uuid:block)) %>% 
  mutate(p = ifelse(is.na(p), 0, p)) %>% 
  filter(!is.na(ratio)) %>% 
  pivot_wider(
    names_from = ratio, names_prefix = "r_", 
    values_from = p, 
    id_cols = c(uuid, condition, block, interval)
  ) %>% 
  mutate(f_interval = ifelse(interval == "0", "NO", 
                             ifelse(interval == "8" | interval == "9", "H", 
                                    interval))
  ) %>% 
  mutate(f_interval = factor(f_interval))
sims <- sims %>% ungroup

model_b <- abcrf(
  formula = f_interval ~ r_1 + r_2 + r_3 + r_4 + r_5 + r_6 + r_7, 
  data=sims %>% filter(condition == "binary"),
  ntree = 4e3, 
  max.depth = NULL,
  paral=T
)
model_r <- abcrf(
  formula = f_interval ~ r_1 + r_2 + r_3 + r_4 + r_5 + r_6 + r_7, 
  data=sims %>% filter(condition == "range"),
  ntree = 4e3, 
  max.depth = NULL,
  paral=T
)
model_n <- abcrf(
  formula = f_interval ~ r_1 + r_2 + r_3 + r_4 + r_5 + r_6 + r_7, 
  data=sims %>% filter(condition == "naturalistic"),
  ntree = 4e3, 
  max.depth = NULL,
  paral=T
)

recovery_results <- function(model){
  CM <- model$model.rf$confusion.matrix
  lvls <- colnames(CM)[1:(length(colnames(CM)) - 1)]
  model$model.rf$confusion.matrix %>% 
    as_tibble() %>% 
    select(-class.error) %>% 
    mutate(true = factor(lvls, levels=lvls)) %>% 
    pivot_longer(c(everything(), -true), names_to = "predicted", values_to = "count") %>% 
    mutate(predicted=factor(predicted, levels=lvls))
}

results <- rbind(
  recovery_results(model_b) %>% mutate(domain = "binary"),
  recovery_results(model_r) %>% mutate(domain = "range"),
  recovery_results(model_n) %>% mutate(domain = "naturalistic")
)
accuracy <- function(results, group=NULL){
  results %>% 
    group_by(correct=true == predicted, {{group}}) %>% 
    summarise(count = sum(count)) %>% 
    group_by({{group}}) %>% 
    mutate(p = count / sum(count)) %>% 
    filter(correct)
}
accuracy(results, true)
accuracy(results, domain)

save(model_b, model_r, model_n, results, sims,
     file = "output/ratio-measures/recovery/abcrf_models.RData")


# Explore recovery results ------------------------------------------------
load(file = "output/ratio-measures/recovery/abcrf_models.RData")
# Accuracy when NO
results %>% 
  mutate(correct = true==predicted) %>% 
  filter(true == "NO") %>% 
  group_by(domain) %>% 
  mutate(N = sum(count)) %>% 
  filter(correct) %>% 
  mutate(p = count / N) %>% 
  ungroup %>% 
  summarise(sum(count) / sum(N)) %>% 
  knitr::kable()
# Accuracy when yes
results %>% 
  mutate(correct = true==predicted) %>% 
  filter(true != "NO") %>% 
  group_by(true, domain) %>% 
  mutate(N = sum(count)) %>% 
  filter(correct) %>% 
  mutate(p = count / N) %>% 
  ungroup %>% 
  summarise(sum(count) / sum(N)) %>% 
  knitr::kable()
# Accuracy when yes if +- 1
results %>% 
  filter(true %in% c(2:6)) %>% 
  mutate(across(c(true, predicted), \(x){as.numeric(as.character(x))})) %>% 
  mutate(correct = ifelse(is.na(predicted), F, abs(true - predicted) <= 1)) %>% 
  group_by(domain, true) %>% 
  mutate(N = sum(count)) %>% 
  filter(correct) %>% 
  summarise(count=sum(count), N = mean(N)) %>% 
  mutate(p = count / N) %>% 
  ungroup %>% 
  summarise(sum(count) / sum(N)) %>% 
  knitr::kable()
# accuracy per domain
accuracy(results, domain)
# conservative / lenient when wrong?
results %>% 
  mutate(correct = true == predicted) %>% 
  filter(!correct) %>% 
  mutate(
    conservative = true != "NO" & (as.numeric(true) < as.numeric(predicted) | predicted == "NO")
  ) %>% 
  group_by(conservative) %>% 
  summarise(count = sum(count)) %>% 
  mutate(p = count / sum(count))
