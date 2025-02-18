rm(list = ls())
library(abcrf)
library(colorspace)
library(glue)
library(lmerTest)
library(rstanarm)
library(patchwork)
library(merTools)
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores() - 1)

select <- dplyr::select
options(mc.cores = parallel::detectCores())
set.seed(2024)
options(contrasts = c("contr.helmert", "contr.poly"))
theme_set(theme_minimal(18))


measure_names <- c("R", "A", "TPF", "D", "RNG", "SWLZ_ER")
threshold <- qnorm(.95)
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
load("output/ratio-measures/recovery/abcrf_models.RData")
data <- get(load("output/ratio-measures/people_p.RData")) %>% 
  filter(exp %in% exps)


# Add domain, study and output to df
data <- data %>% 
  nest(.by=exp) %>% 
  mutate(domain = ifelse(str_detect(exp, "guseva|groot"), "binary", 
                         ifelse(str_detect(exp, "cooper|towse|castillo24_e2"), "range",
                                "naturalistic"))
  ) %>% 
  mutate(study = sapply(exp, \(x){str_split_1(x, "_")[1]})) %>% 
  mutate(output = ifelse(str_detect(study, "groot|guseva"), "keypress",
                         ifelse(str_detect(study, "cooper"), "click", 
                                "voice"))
  ) %>% 
  unnest(data)
# substitute EP=NA for EP=0, 
# remove full sequences (ratio=NA), 
# and pivot wider
data2 <- data %>% 
  mutate(p = ifelse(is.na(p), 0, p)) %>% 
  filter(!is.na(ratio)) %>% 
  pivot_wider(
    names_from = ratio, names_prefix = "r_", 
    values_from = p, 
    id_cols = c(uuid, condition, block, domain, exp, bpm, study, output)
  ) %>% 
  ungroup
data2 <- data2 %>% 
  arrange(uuid, condition, block)


# Get training history ----------------------------------------------------
source("src/rforest_memories.R")
mem_b <- readRDS("output/ratio-measures/recovery/mem_b.rds")
mem_r <- readRDS("output/ratio-measures/recovery/mem_r.rds")
mem_n <- readRDS("output/ratio-measures/recovery/mem_n.rds")




## get predictions (allocations) from forests
pred_b <- predict(
  model_b, 
  obs = data2 %>% filter(domain == "binary"), 
  training = sims %>% filter(condition == "binary"),
  ntree = model_b$model.rf$num.trees, 
  paral = T
)
pred_r <- predict(
  model_r, 
  obs = data2 %>% filter(domain == "range"), 
  training = sims %>% filter(condition == "range"),
  ntree = model_r$model.rf$num.trees, 
  paral = T
)
pred_n <- predict(
  model_n, 
  obs = data2 %>% filter(domain == "naturalistic"), 
  training = sims %>% filter(condition == "naturalistic"), 
  ntree = model_n$model.rf$num.trees, 
  paral = T
)

opts <- levels(pred_n$allocation)

get_posterior <- function(dat, mem, abcrf_model){
  dat.p <- cbind(dat, predict(abcrf_model$model.lda, dat)$x)
  compute_posterior <- function(row){
    leaves <- sapply(1:abcrf_model$model.rf$num.trees, FUN = \(k){find_leaf(treeInfo(abcrf_model$model.rf, k), datum = dat.p[row,])[1]})
    RES <- t(sapply(1:abcrf_model$model.rf$num.trees, FUN = \(x){mem[[x]][leaves[x],]}))
    RES.n <- t(apply(RES, 1, \(x){x/sum(x)}))
    apply(RES.n, 2, sum) / nrow(RES.n)
  }
  
  temp <- foreach(
    r=1:nrow(dat.p), 
    .combine = "rbind", 
    .packages = c("ranger", "abcrf"),
    .export = c("find_leaf")
  ) %dopar% {
    compute_posterior(r)
  }
  
  temp %>% 
    magrittr::set_rownames(NULL) %>% 
    magrittr::set_colnames(str_c("M_", opts)) %>% 
    as_tibble
}


# Get posterior (set to T to do, or F to skip and load, as very time consuming) 
if (F){
  posterior_b <- get_posterior(
    dat=data2 %>% filter(domain=="binary"),
    mem=mem_b,
    abcrf_model = model_b
  )
  saveRDS(posterior_b, "output/ratio-measures/posterior_b.rds")
  posterior_r <- get_posterior(
    dat=data2 %>% filter(domain=="range"),
    mem=mem_r,
    abcrf_model = model_r
  )
  saveRDS(posterior_r, "output/ratio-measures/posterior_r.rds")
  posterior_n <- get_posterior(
    dat=data2 %>% filter(domain=="naturalistic"),
    mem=mem_n,
    abcrf_model = model_n
  )
  saveRDS(posterior_n, "output/ratio-measures/posterior_n.rds")
} else{
  posterior_b <- readRDS("output/ratio-measures/posterior_b.rds")
  posterior_r <- readRDS("output/ratio-measures/posterior_r.rds")
  posterior_n <- readRDS("output/ratio-measures/posterior_n.rds")
}



all_posteriors <- bind_rows(
  data2 %>% 
    filter(domain == "binary") %>% 
    bind_cols(posterior_b) %>% 
    mutate(allocation = pred_b$allocation),
  data2 %>% 
    filter(domain == "range") %>% 
    bind_cols(posterior_r) %>% 
    mutate(allocation = pred_r$allocation),
  data2 %>% 
    filter(domain == "naturalistic") %>% 
    bind_cols(posterior_n) %>% 
    mutate(allocation =pred_n$allocation)
)


# deterministic strategies allocations
det_strat_posteriors <- all_posteriors %>% 
  filter(str_detect(exp, "det_")) %>% 
  select(exp, uuid, allocation) %>% 
  separate(uuid, sep="\\.", into = c("det", "domain", "strategy", "NN", "sim")) %>% 
  mutate(across())

det_strat_posteriors %>% 
  group_by(exp, allocation) %>% 
  tally %>% 
  mutate(p = n / sum(n)) %>% 
  mutate(round(p * 100, 1)) %>% 
  knitr::kable()

det_strat_posteriors %>% 
  group_by(exp, allocation) %>% 
  tally %>% 
  group_by(allocation) %>% 
  summarise(n = sum(n)) %>% 
  mutate(p = n / sum(n) * 100)

#   |exp              |allocation |    n|        p| round(p * 100, 1)|
#   |:----------------|:----------|----:|--------:|-----------------:|
#   |det_binary       |NO         | 7100| 0.887500|              88.8|
#   |det_binary       |H          |  900| 0.112500|              11.2|
#   |det_naturalistic |NO         | 7119| 0.889875|              89.0|
#   |det_naturalistic |6          |   23| 0.002875|               0.3|
#   |det_naturalistic |7          |    1| 0.000125|               0.0|
#   |det_naturalistic |H          |  857| 0.107125|              10.7|
#   |det_range        |NO         | 7723| 0.965375|              96.5|
#   |det_range        |H          |  277| 0.034625|               3.5|

# Percentage of participants who become random at some point -----------------
all_posteriors %>% 
  group_by(study, exp, condition, allocation == "NO") %>% 
  tally %>% 
  mutate(pct = n / sum(n)) %>% 
  filter(!`allocation == "NO"`) %>% 
  filter(!(exp %in% c("half_speed", "det_binary", "det_range", "det_naturalistic"))) %>% 
  ungroup %>% 
  mutate(N = n / pct) %>% 
  mutate(global_pct = sum(n * pct) / sum(N)) %>% 
  mutate(global_N = sum(n*pct))

# Histogram of allocations --
all_posteriors %>% 
  filter(!(exp %in% c("half_speed", "det_binary", "det_range", "det_naturalistic"))) %>% 
  group_by(allocation, domain) %>% 
  tally %>% 
  group_by(domain) %>% 
  mutate(n = n / sum(n)) %>% 
  ggplot(aes(allocation, n, fill=domain)) +
  geom_col(position="dodge")


# pct allocation=NO ~ domain + output -------------------------------------
pct_all_no <- all_posteriors %>% 
  filter(!(exp %in% c("half_speed", "det_binary", "det_range", "det_naturalistic"))) %>% 
  mutate(never_random = allocation == "NO") %>% 
  group_by(exp, output, domain) %>% 
  summarise(never_random = mean(never_random) + 1e-8) %>% 
  mutate(nr.l = boot::logit(never_random))

pct_all_no <- pct_all_no %>% 
  ungroup %>% 
  mutate(across(c(output, domain), factor))

model <- lm(
    data=pct_all_no, formula = nr.l ~ output * domain
  )
anova(model)  

BayesFactor::anovaBF(formula = nr.l ~ output + domain, data=pct_all_no)


## df with estimated allocation (weighted sum of posterior)
times <- all_posteriors %>% 
  filter(allocation != "NO") %>% 
  pivot_longer(M_NO:M_H) %>% 
  filter(name != "M_NO") %>% 
  nest(.by=name) %>% 
  rowwise() %>% 
  mutate(r = str_split_1(name, "_")[2]) %>% 
  mutate(r = ifelse(r == "H", "9", r)) %>% 
  unnest(data) %>% 
  group_by(across(uuid:output)) %>% 
  mutate(valuep = value / sum(value)) %>% 
  mutate(rp = as.numeric(r) * valuep) %>% 
  summarise(rp = sum(rp)) %>% 
  mutate(spi = 60/bpm * rp) %>% 
  ungroup
times2 <- times
times <- times %>% filter(!(exp %in% c("half_speed", "det_binary", "det_range", "det_naturalistic")))

# Summary per exp. condition
pooled_df <- times %>% 
  # filter(study != "speed") %>% 
  group_by(study, exp, condition) %>% 
  mutate(BPM = mean(bpm)) %>% 
  group_by(study, exp, domain, output, condition, BPM) %>% 
  do(
    {model <- lm(spi ~ 1, data = .)
    summary(model)$coefficients %>% 
      as_tibble() %>% 
      magrittr::set_colnames(c("Est", "SE", "t", "p"))}
  ) %>% 
  arrange(Est)
saveRDS(pooled_df, "data/plot_data/spi_pooled_df.rds") # save for plot
# quick plot of spread
pooled_df %>% 
  ungroup %>% 
  mutate(bpm = cut_number(BPM, 3)) %>% 
  mutate(lb = Est - qnorm(0.975) * SE) %>% 
  mutate(ub = Est + qnorm(0.975) * SE) %>% 
  ggplot(aes(forcats::fct_reorder(condition, -Est), Est)) + 
  geom_pointrange(aes(ymin=lb, ymax=ub, fill=domain, color=domain, shape=output, size=bpm)) + 
  geom_hline(yintercept = mean(pooled_df$Est), linetype="dashed") + 
  scale_size_manual(values=1:3 / 3) + 
  coord_flip()

mean(pooled_df$Est)
## Model with output and domain
times <- times %>% mutate(across(c(output, domain, uuid), factor))
model <- lmer(
  data=times, formula = spi ~ output + domain + (1|uuid)
)
if (F){ # time consuming - set to T to carry out
  model.b <- rstanarm::stan_glmer(
    data=times, formula = spi ~ output + domain + (1|uuid)
  )
  bayestestR::bayesfactor_parameters(model.b)
  # Parameter   |       BF
  # ----------------------
  #   (Intercept) | 1.25e+22
  # output1     |    0.656
  # output2     | 4.11e+04
  # domain1     |    0.033
}
set.seed(2024)  # I set seed again here so that results below are same irrespective of 
                # whether the bayesian model ran or was skipped

contrasts(times$output)
summary(model)
car::Anova(model)

# get emmeans for each output type
model2 <- model
model2@beta[2] <- 0 # non-credible
model2@beta[4] <- 0 # non-credible

times %>% 
  distinct(uuid, output, domain) %>% 
  mutate(pred = predict(model2, .)) %>% 
  group_by(output) %>% 
  summarise(mean(pred))


# Variability among participants ------------------------------------------
random_intercepts <- coefficients(model)$uuid[,1]
IQR(random_intercepts)

times %>% 
  select(spi, output, domain, uuid) %>% 
  mutate(domain = factor(domain, levels=c("binary", "range", "naturalistic"))) %>%
  mutate(pred = predict(model, .)) %>% 
  mutate(IQR = IQR(pred)) %>% 
  mutate(lb = pred - .5 * IQR, ub = pred + .5 * IQR) %>% 
  summarise(across(pred:ub, mean))

# speed -------------------------------------------------------------------
plot_speed <- function(var){
  speed_times %>% 
    mutate(condition = str_remove(condition, "speed.")) %>% 
    pivot_wider(names_from = condition, values_from = {{var}}, id_cols = c(uuid, fastFirst)) %>% 
    drop_na() %>% 
    ggplot(aes(Fast, Slow)) + 
    geom_point(aes(fill=fastFirst), size=3, alpha=.5, shape=21, color="black") +
    geom_function(fun=\(x){x}) + 
    geom_function(fun=\(x){1/2*x}, linetype="dashed")
}
identity_error <- function(var){
  temp <- speed_times %>% 
    mutate(condition = str_remove(condition, "speed.")) %>% 
    pivot_wider(names_from = condition, values_from = {{var}}, id_cols = c(uuid, fastFirst)) %>% 
    drop_na()
  SSE <- temp %>% 
    mutate(error = Slow - mean(Slow)) %>% 
    mutate(sqerror = error**2) %>% 
    summarise(SSE=sum(sqerror)) %>% 
    pull(SSE)
  SSE2 <- temp %>% 
    mutate(error = Slow - Fast) %>% 
    mutate(sqerror = error**2) %>% 
    summarise(SSE=sum(sqerror)) %>% 
    pull(SSE)
  
  SSE2
}
predict_using_rp <- function(
    bpm_Fast, 
    bpm_Slow, 
    rp_Fast
){
  rp_Fast
}
predict_using_spi <- function(
    bpm_Fast, 
    bpm_Slow, 
    rp_Fast
){
  spi_Fast <- 60 / bpm_Fast * rp_Fast
  spi_Fast * bpm_Slow / 60
}
filter_speed_times <- function(speed_exp){
  times2 %>% 
    filter(str_detect(exp, {{speed_exp}})) %>% 
    mutate(condition = ifelse(str_detect(condition, "Fast"), "Fast", "Slow")) %>% 
    mutate(block = factor(block), condition = factor(condition)) %>% 
    mutate(fastFirst = (str_detect(condition, "Fast") & block == 1) | (str_detect(condition, "Slow") & block==2))
}

speed_exp <- c("speed2", "half_speed", "towse16")[1] # change this to do with TOWSE (Appendix B)

# filter times to only this df + add fastFirst condition
speed_times <- filter_speed_times(speed_exp)

# Save dfs for possible plotting
for (speed_exp in c("speed2", "half_speed", "towse16")){
  saveRDS(filter_speed_times(speed_exp), glue("data/plot_data/speedtimes_{speed_exp}.rds"))
}


lik <- function(obs, exp, mean, sd){
  den <- sqrt(2 * mean * sd**2)
  expnt <- .5 * ((obs-exp) / sd)**2
  (1 / den) * exp(-expnt)
}

speed_times %>%
  pivot_wider(names_from = condition, values_from = c(bpm, rp, spi), id_cols = uuid) %>%
  drop_na() %>% 
  select(uuid, rp_Fast, rp_Slow) %>% 
  mutate(mod_1 = rp_Fast) %>% 
  mutate(mod_half = 1/2 * rp_Fast) %>% 
  mutate(M = mean(rp_Slow)) %>% 
  mutate(SD = sd(rp_Slow)) %>% 
  mutate(L_mod_1 = lik(rp_Slow, mod_1, M, SD)) %>% 
  mutate(L_mod_half = lik(rp_Slow, mod_half, M, SD)) %>% 
  pivot_longer(L_mod_1:L_mod_half, names_to = "model") %>% 
  group_by(model) %>% 
  summarise(AIC = -2 * sum(log(value))) %>% 
  mutate(AICw = cogmod::AICw(AIC)[[2]])

saveRDS(speed_times, glue("data/plot_data/speedtimes_{speed_exp}.rds")) # save for data
