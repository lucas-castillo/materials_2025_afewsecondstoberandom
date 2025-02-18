rm(list=ls())
library(tidyverse)
library(glue)
library(lmerTest)
library(rstanarm)
library(logspline)
library(bayestestR)
options(mc.cores=parallel::detectCores() - 1)
set.seed(2023)
odds <- function(p){return(p / (1-p))}
min_iter <- 6000

for (exp in c("fs")){
  # Setup ---
  dir.create(
    glue("output/descriptive-analyses/{exp}/"), 
    showWarnings = F, recursive = T
  )
  data <- get(load(glue("data/rdata/{exp}_fulldescriptives.RData"))) 
  data_sv <- get(load(glue("data/rdata/{exp}_sv.RData"))) 
  if (exp == "fs"){
    data.t <- get(load(glue("data/rdata/fs_fulldescriptives_trunc.RData"))) 
    data_sv.t <- get(load(glue("data/rdata/fs_sv_trunc.RData"))) 
    contrasts(data.t$condition) <- c(-1, 1)
    contrasts(data_sv.t$condition) <- c(-1, 1)
    
    contrasts(data.t$blockorder) <- c(-1, 1)
    contrasts(data_sv.t$blockorder) <- c(-1, 1)
    
    contrasts(data.t$block) <- c(-1, 1)
    contrasts(data_sv.t$block) <- c(-1, 1)
  }
  
  contrasts(data$condition) <- c(-1, 1)
  contrasts(data_sv$condition) <- c(-1, 1)
  
  contrasts(data$blockorder) <- c(-1, 1)
  contrasts(data_sv$blockorder) <- c(-1, 1)
  
  contrasts(data$block) <- c(-1, 1)
  contrasts(data_sv$block) <- c(-1, 1)
  
  measure_names <- data %>% 
    select(contains("boot_")) %>% 
    colnames %>% 
    str_remove("boot_")  
  logit_measures <- c("R", "A", "TPF")
  sv_measures <- data_sv %>% 
    select(-c(uuid:blockorder), -contains("boot")) %>% 
    colnames
  ## Do expectations depend on condition? ---------------
  expectations <- data %>% 
    distinct(uuid, condition, .keep_all = T)
  
  saved_object_names <- c()
  for (measure in measure_names){
    # Make formula for models
    if (measure %in% logit_measures){
      eval(parse(text=glue("formula_{measure} <- I(boot::logit(boot_{measure})) ~ condition + (1|uuid)")))
    } else {
      eval(parse(text=glue("formula_{measure} <- boot_{measure} ~ condition + (1|uuid)")))
    }
    # make frequentist model
    eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = expectations)")))
    # Make bayesian model
    eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = expectations, iter = min_iter)")))
    # get bf
    eval(parse(text=glue("mod{measure}.bf <- bayesfactor_parameters(mod{measure}.b, null = 0)")))
    # add objects to list to save
    saved_object_names <- c(
      saved_object_names, 
      glue("mod{measure}"),
      glue("mod{measure}.b"),
      glue("mod{measure}.bf")
    )
  }
  # save to file  
  save(
    list=saved_object_names,
    file = glue("output/descriptive-analyses/{exp}/expectations-condition.RData")
  )

  ## Are people different to expectations? (i.e. vs Intercept) -----
  saved_object_names <- c()
  for (measure in measure_names){
    # Make formula for models
    if (measure %in% logit_measures){
      eval(parse(text=glue("formula_{measure} <- {measure} ~ offset(log(odds(boot_{measure}))) + (1|uuid)")))
    } else {
      eval(parse(text=glue("formula_{measure} <- I({measure} - boot_{measure}) ~ (1|uuid)")))
    }
    # make frequentist model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure} <- glmer(formula_{measure}, data = data, family='binomial')")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data_sv)")))
    } else {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data)")))
    }
    # Make bayesian model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, family='binomial', iter=min_iter)")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data_sv, iter=min_iter)")))
    } else {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, iter=min_iter)")))
    }
    # get bf
    eval(parse(text=glue("mod{measure}.bf <- bayesfactor_parameters(mod{measure}.b, null = 0)")))
    # add objects to list to save
    saved_object_names <- c(
      saved_object_names, 
      glue("mod{measure}"),
      glue("mod{measure}.b"),
      glue("mod{measure}.bf")
    )
  }
  # save to file  
  save(
    list=saved_object_names,
    file = glue("output/descriptive-analyses/{exp}/intercept.RData")
  )

  ## Measures vs condition (no need for bootstrap offset) -----
  saved_object_names <- c()
  for (measure in measure_names){
    # Make formula for models
    eval(parse(text=glue("formula_{measure} <- {measure} ~ condition + (1|uuid)")))
    
    data_names <- c("data", "data_sv")
    if (exp == "fs") data_names <- str_c(data_names, ".t")
    
    # make frequentist model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure} <- glmer(formula_{measure}, data = {data_names[1]}, family='binomial')")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = {data_names[2]})")))
    } else {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = {data_names[1]})")))
    }
    # Make bayesian model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = {data_names[1]}, family='binomial', iter=min_iter)")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = {data_names[2]}, iter=min_iter)")))
    } else {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = {data_names[1]}, iter=min_iter)")))
    }
    # get bf
    eval(parse(text=glue("mod{measure}.bf <- bayesfactor_parameters(mod{measure}.b, null = 0)")))
    # add objects to list to save
    saved_object_names <- c(
      saved_object_names, 
      glue("mod{measure}"),
      glue("mod{measure}.b"),
      glue("mod{measure}.bf")
    )
  }
  # save to file  
  save(
    list=saved_object_names,
    file = glue("output/descriptive-analyses/{exp}/condition.RData")
  )
  
  # Order x Condition --------------------
  saved_object_names <- c()
  for (measure in measure_names){
    # Make formula for models
    eval(parse(text=glue("formula_{measure} <- {measure} ~ blockorder * condition + (1|uuid)")))
    
    # make frequentist model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure} <- glmer(formula_{measure}, data = data, family='binomial')")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data_sv)")))
    } else {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data)")))
    }
    # Make bayesian model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, family='binomial', iter=min_iter)")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data_sv, iter=min_iter)")))
    } else {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, iter=min_iter)")))
    }
    # get bf
    eval(parse(text=glue("mod{measure}.bf <- bayesfactor_parameters(mod{measure}.b, null = 0)")))
    # add objects to list to save
    saved_object_names <- c(
      saved_object_names, 
      glue("mod{measure}"),
      glue("mod{measure}.b"),
      glue("mod{measure}.bf")
    )
  }
  # save to file  
  save(
    list=saved_object_names,
    file = glue("output/descriptive-analyses/{exp}/order_condition.RData")
  )
  
  # Block x Condition --------------------
  saved_object_names <- c()
  for (measure in measure_names){
    # Make formula for models
    eval(parse(text=glue("formula_{measure} <- {measure} ~ block * condition + (1|uuid)")))
    
    # make frequentist model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure} <- glmer(formula_{measure}, data = data, family='binomial')")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data_sv)")))
    } else {
      eval(parse(text=glue("mod{measure} <- lmer(formula_{measure}, data = data)")))
    }
    # Make bayesian model
    if (measure %in% logit_measures){
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, family='binomial', iter=min_iter)")))
    } else if (measure %in% sv_measures) {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data_sv, iter=min_iter)")))
    } else {
      eval(parse(text=glue("mod{measure}.b <- stan_glmer(formula_{measure}, data = data, iter=min_iter)")))
    }
    # get bf
    eval(parse(text=glue("mod{measure}.bf <- bayesfactor_parameters(mod{measure}.b, null = 0)")))
    # add objects to list to save
    saved_object_names <- c(
      saved_object_names, 
      glue("mod{measure}"),
      glue("mod{measure}.b"),
      glue("mod{measure}.bf")
    )
  }
  # save to file  
  save(
    list=saved_object_names,
    file = glue("output/descriptive-analyses/{exp}/block_condition.RData")
  )
}
