library(tidyverse)
library(glue)
source("src/get_1D_summaries.R")

boxcox_findlambda <- function(v){
  bc <- MASS::boxcox(v ~ 1, plotit=F)
  lambda <- bc$x[which(bc$y == max(bc$y))]
  return(lambda[1])
}

boxcox_transform <- function(v, lambda){
  if (lambda != 0){
    v <- (v ** (lambda) - 1) / lambda
  } else(
    v <- log(v)
  )
  return(v)
}
non_overlap_seqs <- function(
    N=NULL, 
    max_ratio=7,
    equal_numbers=T
  ){
  if (is.null(N)){
    stop("Please provide a sequence length")
  }
  sequence_length <- floor(N / max_ratio)
  legal_starts <- tibble()
  for (rat in 1:max_ratio){
    legal_indices <- matrix(ncol = sequence_length)
    continue <- T
    k <- 0
    while (continue){
      add <- sapply(
        k * sequence_length * rat + 1:rat, ## starting indices
        \(i) {seq(i, length.out=sequence_length, by=rat)}) %>% 
        t()
      k <- k + 1
      for (r in 1:nrow(add)){
        if (all(add[r,] < N)){
          legal_indices <- rbind(legal_indices, add[r,])
        }
        else{
          continue <- F    
        }
      }
    }
    legal_starts <- rbind(
      legal_starts, 
      tibble(
        i=legal_indices[!is.na(legal_indices[,1]),1], 
        ratio=rat, 
        N=N, 
        sequence_length
      ))
  }
  if (!equal_numbers){
    return(legal_starts)
  }
  return(
    legal_starts %>% 
      group_by(ratio) %>% 
      mutate(min_ratio = n(), idx = 1:n()) %>% 
      ungroup() %>% 
      mutate(min_ratio = min(min_ratio)) %>% 
      filter(idx <= min_ratio) %>% 
      dplyr::select(-c(min_ratio, idx))
  )
}
compute_ratio_thinning <- function(
  v,
  ratios,
  do_iid=T,
  N_IID=2e3, 
  is.circular=F, 
  drop.r=F
){
  # fix sequence if needed
  v <- v[!is.na(v)]
  v <- round(v)
  # get n alternatives
  lb <- min(v)
  ub <- max(v)
  alt <- lb:ub
  # get length and thinned length
  n <- length(v)
  thinned_length <- floor(n / max(ratios))
  
  sim <- data.frame()
  # get measures for the sequence
  sim <- rbind(
    sim,
    all_measures(
      v, 
      alt, 
      drop.r, 
      is.circular=is.circular
    ) %>% 
      mutate(type="full", source="participant", N=n, ratio=NA)
  )
  # For each ratio, which starting i make independent sequences (no shared items)
  # so that each ratio has the same number of starting indices
  start_indices <- non_overlap_seqs(N=n, max_ratio = max(ratios))
  # get measures for different ratios -------------
  for (r in ratios){
    is <- start_indices %>% 
      filter(ratio == r) %>% 
      pull(i)
    for (i in is){
      indices <- seq(from=i, by=r, length.out = thinned_length)
      temp_v <- v[indices]
      sim <- rbind(
        sim,
        all_measures(
          temp_v, 
          alt, 
          drop.r, 
          is.circular=is.circular
        ) %>% 
          mutate(
            type={if (r == 1) glue("wind_{i}") else glue("thin_{i}")}, 
            source="participant", 
            N=thinned_length, 
            ratio=r
          )
      )
    }
  }
  
  # get iid simulations for N=n and N=thinned_length
  if (do_iid){
    for (i in 1:N_IID){
      sim <- rbind(
        sim, 
        all_measures(
          sample(v), 
          alt, 
          drop.r, 
          is.circular=is.circular
        ) %>% mutate(
          type="full", 
          source=glue("iid_{str_pad(i, nchar(N_IID), pad = 0)}"), 
          N=n, 
          ratio=NA
          ),
        all_measures(
          sample(v)[1:thinned_length], 
          alt, 
          drop.r, 
          is.circular=is.circular
        ) %>% 
          mutate(
            type="thin", 
            source=glue("iid_{str_pad(i, nchar(N_IID), pad = 0)}"), 
            N=thinned_length, 
            ratio=NA
          )
      )
    }
  }
  return(sim)
}

process_ratio_measures <- function(
    M,
    measure_names,
    uuid,
    condition,
    block,
    exp
){
  if (exp == "cooper12"){
    measure_names <- measure_names[measure_names != "TPF"]
  }
  temp <- M %>%
    mutate(sequence_id = str_c(uuid, condition, block, sep = "__")) %>% 
    mutate(uuid = uuid, condition = condition, block = block)
  
  
  if (exp == "cooper12"){
    temp <- temp %>% 
      select(-TPF)
  }
  ## Calc norm value for each measure -----------------------------------------
  temp <- temp %>% 
    ungroup %>%
    pivot_longer(any_of(measure_names), names_to="measure") %>% 
    group_by(measure, N) %>% 
    mutate(lambda = boxcox_findlambda(value + 1e-8)) %>% 
    rowwise() %>% 
    mutate(norm_value = boxcox_transform(value + 1e-8, lambda)) %>% 
    group_by(measure, N) %>% 
    mutate(sd_norm = sd(norm_value, na.rm = T)) %>% 
    mutate(norm_value = (norm_value - mean(norm_value, na.rm=T)) / sd_norm) %>% 
    ungroup()
  
  
  # Safety: remove measures with sd = 0 -------------------------------------
  bad_measures <- temp %>% 
    group_by(measure, N) %>% 
    summarise(sd_norm = mean(sd_norm)) %>% 
    filter(sd_norm == 0) %>% 
    pull(measure)
  
  bad_measures <- c(
    bad_measures, 
    temp %>% 
      group_by(source == "participant", measure) %>% 
      summarise(problem=all(is.na(norm_value))) %>% 
      filter(problem) %>% 
      pull(measure)
  )
  
  # measure_names <- measure_names[!(measure_names %in% bad_measures)]
  
  temp <- temp %>% 
    filter(measure %in% measure_names)
  # split data into people and iid ------------
  temp_person <- temp %>%
    filter(source == "participant")
  temp_iid <- temp %>%
    filter(source != "participant")
  # Store mean + sd of iid values per uuid, condition, measure and N (thresholds) --------
  temp_iid_summaries <- temp_iid %>%
    pivot_wider(
      names_from = measure, values_from = norm_value,
      id_cols = c(uuid, condition, ratio, type, source, N)
    ) %>%
    group_by(N) %>%
    summarise_at(.vars = measure_names,
                 .funs = list("M"=mean, "SD"=sd), na.rm=T) %>%
    pivot_longer(-c(N), names_to = "measure") %>%
    mutate(summary = ifelse(str_detect(measure, "_M"), "M", "SD")) %>%
    mutate(measure = str_remove(measure, "_M")) %>%
    mutate(measure = str_remove(measure, "_SD")) %>%
    pivot_wider(names_from = summary, values_from = value) %>%
    ungroup
  
  # Transfer iid summaries to people + iid dataframe ------------------------------
  temp_person <- temp_person %>%
    rowwise() %>%
    mutate(ivalue = temp_iid_summaries$M[
      temp_iid_summaries$N == N &
        temp_iid_summaries$measure == measure
    ]) %>%
    mutate(isd = temp_iid_summaries$SD[
      temp_iid_summaries$N == N &
        temp_iid_summaries$measure == measure
    ])
  
  temp_iid <- temp_iid %>%
    nest(.by = c(N, measure)) %>%
    rowwise() %>%
    mutate(ivalue = temp_iid_summaries$M[
      temp_iid_summaries$N == N &
        temp_iid_summaries$measure == measure
    ]) %>%
    mutate(isd = temp_iid_summaries$SD[
      temp_iid_summaries$N == N &
        temp_iid_summaries$measure == measure
    ]) %>%
    unnest(data) %>%
    ungroup()
  
  # Add distance measure d to people + iid dfs ----------------------------
  temp_person <- temp_person %>%
    ungroup %>%
    mutate(d = norm_value - ivalue)
  
  temp_iid <- temp_iid %>%
    ungroup %>%
    mutate(d = norm_value - ivalue)
  # Separate source and simulation number  -------------------------------------
  temp_person <- temp_person %>%
    separate(type, into=c("type", "i"), sep="_", fill="right")
  temp_iid <- temp_iid %>% 
    separate(source, into=c("source", "i"), sep="_", fill="right")
  
  
  # Calculate euclidean distance between each sequence and iid expectation --------
  temp_iid.distance <- temp_iid %>%
    {if (exp == "cooper12") {
      filter(., measure != "TPF") 
    } else {
      filter(.)
    }} %>% 
    group_by(uuid, condition, block, N, type, source, i) %>%
    summarise(
      ed        = sqrt(sum((norm_value - ivalue)**2))
    ) %>%
    ungroup
  
  empirical_p <- function(v, n){
    temp_iid.distance %>% 
      filter(N == n) %>% 
      mutate(surpass = ed >= v) %>% 
      summarise(surpass = mean(surpass)) %>% 
      pull(surpass)
  }
  
  temp_person.p <- temp_person %>% 
    {if (exp == "cooper12") {
      filter(., measure != "TPF") 
    } else {
      filter(.)
    }} %>% 
    group_by(uuid, condition, ratio, type, i, N) %>%
    summarise(
      ed        = sqrt(sum((norm_value - ivalue)**2))
    ) %>% 
    ungroup
  
  temp_person.p <- temp_person.p %>% 
    group_by(uuid, condition, ratio, type) %>% 
    summarise(ed = sum(ed), N=mean(N), n_subseq=n())
  
  n_samples <- temp_person.p %>% 
    filter(!is.na(ratio)) %>% 
    pull(n_subseq) %>% 
    first()
  
  if (length(temp_iid.distance$ed[!is.na(temp_iid.distance$ed)]) == 0){
    bar <- rep(0, 2000)
    bar_1 <- rep(0, 2000)
  } else {
    bar <- sapply(1:2000, FUN = \(x){sum(sample(temp_iid.distance$ed[!is.na(temp_iid.distance$ed)], size = n_samples, replace = T))})
    bar_1 <- sapply(1:2000, FUN = \(x){sum(sample(temp_iid.distance$ed[!is.na(temp_iid.distance$ed)], size = 1))})
    
  }
    
  temp_person.p <- temp_person.p %>% 
    mutate(p = mapply(\(d, n){
      if(n == 1){
        mean(bar_1 >= d)
      }else(
        mean(bar >= d)
      )
    }, d=ed, n=n_subseq))
  
  temp_person.p <- temp_person.p %>% 
    mutate(threshold=sapply(n_subseq, \(x){
      if (x == 1){
        quantile(bar_1, .95, na.rm = T)
      } else{
        quantile(bar, .95, na.rm = T)
      }
    }))
  
  # ensure uuid,condition,block in everything saved ----------------------------
  temp_iid <- temp_iid %>%
    mutate(uuid = uuid, condition = condition, block=block) %>%
    relocate(uuid, condition, block)
  temp_iid.distance <- temp_iid.distance %>%
    mutate(uuid = uuid, condition = condition, block=block) %>%
    relocate(uuid, condition, block)
  temp_iid_summaries <- temp_iid_summaries %>%
    mutate(uuid = uuid, condition = condition, block=block) %>%
    relocate(uuid, condition, block)
  temp_person <- temp_person %>%
    mutate(uuid = uuid, condition = condition, block=block) %>%
    relocate(uuid, condition, block)
  temp_person.p <- temp_person.p %>%
    mutate(uuid = uuid, condition = condition, block=block) %>%
    relocate(uuid, condition, block)
  
  return(list(temp_iid, temp_iid.distance, temp_iid_summaries, temp_person, temp_person.p))
}

get_stochastic_intervals <- function(interval, N, threshold=.25){
  try=T
  while (try){
    interval_values <- if (interval != 1){rep(c(0,1,2), 5) + interval - 1} else c(rep(1, 5), 2)
    planned_intervals <- sample(rep(interval_values, ceiling(N / sum(interval_values))))
    planned_intervals.cumsum <- cumsum(planned_intervals)
    planned_intervals <- planned_intervals[1:min(which(planned_intervals.cumsum >= N))]
    produced_interval <- mean(planned_intervals[1:min(which(planned_intervals.cumsum >= N))])
    if (abs(produced_interval - interval) < threshold){
      try=F
    }
  } 
  return(planned_intervals)
}

generate_random <- function(
    domain = "range",
    det_f = \(x){x + 1},
    N = 100,
    interval = 3
){
  if (interval == 0){
    if (domain == "binary" | domain == "range"){
      k <- if (domain == "binary") 2 else 10
      v <- sample(0:(k-1), 1, replace = T)
      for (i in 2:N){
        v[i] <- det_f(v[i-1])
      }
      v <- v %% k
    } else if (domain == "naturalistic"){
      k <- 220
      v <- round(rnorm(1, 170, 12))
      for (i in 2:N){
        v[i] <- det_f(v[i-1])
      }
      v <- v %% k
      v[v < 100] <- v[v < 100] + 100
    }
  } else{
    intervals <- get_stochastic_intervals(interval, N)
    if (domain == "binary" | domain == "range"){
      k <- if (domain == "binary") 2 else 10
      v <- sample(0:(k-1), N, replace = T)
      M <- matrix(v)
      for (i in 2:max(intervals)){
        M <- cbind(M, sapply(M[, ncol(M)], det_f))
      }
      v <- c()
      for (r in 1:length(intervals)){
        v <- c(v, M[r,1:intervals[r]])
      }
      v <- v %% k
    } else if (domain == "naturalistic"){
      v <- round(rnorm(N, 170, 12))
      M <- matrix(v)
      for (i in 2:max(intervals)){
        M <- cbind(M, sapply(M[, ncol(M)], det_f))
      }
      v <- c()
      for (r in 1:length(intervals)){
        v <- c(v, M[r,1:intervals[r]])
      }
    }
  }
  return(v[1:N])  
}
