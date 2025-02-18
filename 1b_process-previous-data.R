library(tidyverse)
library(glue)

# cooper12 ----------------------------------------------------------------
cooper_files <- list.files("data/cooper12/experiment1/")
cooper_files <- cooper_files[str_detect(cooper_files, "0\\.txt")]
aux_data <- read_table("data/cooper12/experiment1/auxiliary_data.tab") %>% 
  mutate(sex = c("male", "female")[sex])
cooper12 <- tibble()

for (f in cooper_files){
  idn <- as.numeric(str_split(f, "_")[[1]][2])
  gdr <- aux_data$sex[aux_data$partnr == idn]
  cooper12 <- rbind(
    cooper12,
    read_table(
      paste("data/cooper12/experiment1/", f, sep=""), 
      col_names = c("response", "rt")) %>% 
      mutate(id = idn, gender = gdr) %>% 
      mutate(response = LETTERS[response + 1])
  )
}
cooper12 <- cooper12 %>% 
  rename(uuid = "id") %>% 
  mutate(uuid = paste("cooper12", uuid, sep=".")) %>% 
  mutate(condition = "cooper12.a") %>% 
  mutate(value = sapply(response, \(r){which(r == LETTERS)}, USE.NAMES = F))

cooper12 <- cooper12 %>% 
  group_by(uuid) %>% 
  mutate(t = sum(rt), n = max(n())) %>% 
  mutate(bpm = n / t * 60) %>% 
  select(-t) %>% 
  ungroup %>% 
  mutate(block = 1)

save(cooper12, file = "data/rdata/cooper12.RData")

# towse16 -----------------------------------------------------------------
e3_coop <- read_csv(
    file = "data/towse16/PLOS_E3_Cooperation_RawData.csv", 
    col_types = paste(rep("c", 40), collapse="")
    ) %>% 
    pivot_longer(cols = c(everything(), -c(Index, Type, Condition)), names_to = "uuid") %>% 
    filter(Type != "Pair")
e3_negl <- read_csv(
    file = "data/towse16/PLOS_E3_Neglect_RawData.csv", 
    col_types = paste(rep("c", 40), collapse="")) %>% 
    pivot_longer(cols = c(everything(), -c(Index, Type, Condition)), names_to = "uuid") %>% 
  filter(Type != "Pair")

towse16_e3 <- rbind(
  e3_coop, e3_negl
) %>% 
  rename(condition = "Type", index = "Index") %>% 
  select(-Condition) %>% 
  mutate(
    value = as.numeric(value), 
    index = as.numeric(index), 
    # condition = factor(condition),
    # uuid = factor(uuid)
  ) %>% 
  mutate(bpm = ifelse(condition == "Slow", 20, 40)) %>% # pace in bpm
  arrange(uuid, condition, index)

towse16_e3 <- towse16_e3 %>% 
  mutate(uuid = paste("towse16_e3", uuid, sep=".")) %>% 
  mutate(condition = paste("towse16_e3", condition, sep="."))

e1_ids <- read_csv("data/towse16/PLOS_E1.csv") %>% 
  pull(pno) %>% 
  str_remove(".00") %>% 
  tolower()
e2_ids <- read_csv("data/towse16/PLOS_E2.csv") %>% 
  pull(pno) %>% 
  str_remove(".00") %>% 
  tolower()
towse16_e1 <- tibble()
towse16_e2 <- tibble()
for (condition in c("Fast", "Slow")){
  print(condition)
  for (f in list.files(glue("data/towse16/Response_Sequences_E1E2/{condition}/"))){
    id <- str_extract(f, "[0-9]+[a-b]?") # get number + a/b as id
    id <- str_remove(id, "^0+") # remove leading 0
    exp <- if (id %in% e1_ids) "e1" else if (id %in% e2_ids) "e2" else NA
    sequence <- read_table(glue("data/towse16/Response_Sequences_E1E2/{condition}/{f}"), col_names = F) %>% 
      pivot_longer(everything()) %>% 
      pull(value)
    sequence <- if (is.na(sequence[length(sequence)])) sequence[1:(length(sequence)-1)] else sequence
    M <- tibble(
      uuid = id, 
      bpm = if (condition == "Fast") 40 else 20,
      value = sequence, 
      condition = condition
    ) %>% 
      mutate(index = 1:n())
    if (exp == "e1"){
      towse16_e1 <- rbind(towse16_e1, M)
    } else if (exp == "e2"){
      towse16_e2 <- rbind(towse16_e2, M)
    } else{
      warning("I don't know where this participant comes from!")
    }
  }
}
towse16_e1 <- towse16_e1 %>% 
  mutate(uuid = paste("towse16_e1", uuid, sep=".")) %>% 
  mutate(condition = paste("towse16_e1", condition, sep="."))
towse16_e2 <- towse16_e2 %>% 
  mutate(uuid = paste("towse16_e2", uuid, sep=".")) %>% 
  mutate(condition = paste("towse16_e2", condition, sep="."))
## We don't know which sequence people did first
towse16_e1 <- towse16_e1 %>% mutate(block=NA)
towse16_e2 <- towse16_e2 %>% mutate(block=NA)
towse16_e3 <- towse16_e3 %>% mutate(block=NA)

towse16 <- rbind(
  towse16_e1,
  towse16_e2,
  towse16_e3
)
save(towse16_e1, file = "data/rdata/towse16_e1.RData")
save(towse16_e2, file = "data/rdata/towse16_e2.RData")
save(towse16_e3, file = "data/rdata/towse16_e3.RData")

# guseva23 ----------------------------------------------------------------
guseva23 <- read_csv("data/guseva23/dat.csv") %>% 
  mutate(choice = (choice - 1) / -2) %>% # 1 heads 0 tails
  rename(uuid = "subject", value = "choice", index = "trial") %>% 
  mutate(condition = c(
    "Explicit Randomness", 
    "Free Choice", 
    "Irregularity", 
    "Mental Toss", 
    "Perceptual Guessing")[condition]) %>% 
  select(-workerID)

guseva23 <- guseva23 %>% 
  group_by(uuid, block, condition) %>% 
  # Trials are fixation (500ms) - choice - feedback (500ms)
  # so ITI is 500ms in 1st trial 1000ms elsewhere
  mutate(ITI = ifelse(index == 1, 500, 1000)) %>% 
  mutate(rt = (rt + ITI) / 1000) %>% 
  mutate(total_time = sum(rt)) %>% 
  mutate(observed_n = n()) %>% 
  mutate(bpm = observed_n / (total_time / 60)) %>% 
  select(-c(ITI:observed_n)) %>% 
  ungroup


guseva23_er <- guseva23 %>% 
  filter(condition == "Explicit Randomness") %>% 
  rename(task = "condition") %>% 
  mutate(uuid = paste("guseva23_er", uuid, sep=".")) %>% 
  mutate(condition = "guseva23_er.a")
guseva23_irr <- guseva23 %>% 
  filter(condition == "Irregularity") %>% 
  rename(task = "condition") %>% 
  mutate(uuid = paste("guseva23_irr", uuid, sep=".")) %>% 
  mutate(condition = "guseva23_irr.a")
guseva23_mt <- guseva23 %>% 
  filter(condition == "Mental Toss") %>% 
  rename(task = "condition") %>% 
  mutate(uuid = paste("guseva23_mt", uuid, sep=".")) %>% 
  mutate(condition = "guseva23_mt.a")

save(guseva23_er, file = "data/rdata/guseva23_er.RData")
save(guseva23_irr, file = "data/rdata/guseva23_irr.RData")
save(guseva23_mt, file = "data/rdata/guseva23_mt.RData")



# groot22 -----------------------------------------------------------------
groot22 <- tibble()
for (f in list.files("data/groot22/")){
  if (!str_detect(f, "MRI")){
    id <- str_extract(f, "[0-9]+")
    M <- read_csv(glue("data/groot22/{f}"), skip = 1)
    M <- M %>% 
      mutate(new_block = stimulus == "stimulus" & trial == 0) %>% 
      mutate(block = cumsum(new_block)) %>% 
      select(-new_block)
    
    M <- M %>%
      filter(block_type == "random" & stimulus == "tap") %>% 
      mutate(uuid = id) %>%
      rename(value = "response", index = "trial") %>%
      mutate(value = sapply(value, \(v){which(v == unique(M$response)) - 1}, USE.NAMES = F)) %>%
      mutate(index = index + 1) %>%
      select(-block_type, -stimulus)
    
    M <- M %>% 
      group_by(block) %>% 
      mutate(delta_time = time - min(time)) %>% 
      mutate(total_time = max(delta_time)) %>% 
      mutate(observed_n = n()) %>% 
      mutate(bpm = (observed_n - 1) / (total_time) * 60) %>% 
      ungroup() %>% 
      select(-c(delta_time:observed_n))
    
    blocks <- sort(unique(M$block))
    new_block <- sapply(M$block, \(x){which(x == blocks)}, USE.NAMES = F)
    M$block <- new_block
    groot22 <- rbind(groot22, M)
  }
}

groot22 <- groot22 %>% 
  group_by(uuid, block) %>% 
  mutate(time = time - min(time)) %>% 
  ungroup %>% 
  # mutate(across(c(block, uuid), factor)) %>% 
  mutate(condition = "groot22.a")

groot22 <- groot22 %>% 
  mutate(uuid = paste("groot22", uuid, sep = "."))

save(groot22, file = "data/rdata/groot22.RData")

# castillo24_e1 --------------------------------------------------------------
castillo24_e1 <- samplrData::castillo2024.rgmomentum.e1 %>% 
  as_tibble() 


castillo24_e1 <- castillo24_e1 %>% 
  select(id, target_dist, block, unit, value, value_in_units) %>% 
  rename(uuid = "id", condition="target_dist") %>% 
  group_by(uuid, condition) %>% 
  mutate(max_units = length(unique(unit))) %>% 
  # we select only sequences who used 1 unit in the sequence (they didn't alternate between cm and f/in)
  filter(max_units == 1) %>% 
  ungroup %>% 
  select(uuid, condition, block, value_in_units) %>% 
  rename(value = "value_in_units") %>% 
  group_by(uuid, condition) %>% 
  mutate(index = 1:n()) %>% 
  mutate(bpm = n() / 5) %>% 
  ungroup() %>% 
  mutate(block = sapply(block, \(x){which(LETTERS == x)}, USE.NAMES = F))

castillo24_e1 <- castillo24_e1 %>%
  mutate(uuid = paste("castillo24_e1", uuid, sep=".")) %>% 
  mutate(condition = paste("castillo24_e1", condition, sep="."))

save(castillo24_e1, file = "data/rdata/castillo24_e1.RData")

# castillo24_e2--------------------------------------------------------------
castillo24_e2 <- samplrData::castillo2024.rgmomentum.e2 %>% 
  as_tibble() 


castillo24_e2 <- castillo24_e2 %>% 
  filter(dim == 1) %>% 
  mutate(block = sapply(block, \(x){which(LETTERS == x)}, USE.NAMES = F)) %>% 
  rename(uuid = "id", condition="dim", value="position") %>% 
  select(uuid, condition, block, value) %>% 
  group_by(uuid, condition, block) %>% 
  mutate(index = 1:n()) %>% 
  mutate(bpm = n() / 5) %>% 
  ungroup

castillo24_e2 <- castillo24_e2 %>%
  mutate(uuid = paste("castillo24_e2", uuid, sep=".")) %>% 
  mutate(condition = paste("castillo24_e2", condition, sep="."))
save(castillo24_e2, file = "data/rdata/castillo24_e2.RData")


# speed -------------------------------------------------------------------
# this only adds bpm to speed
speed <- get(load("data/rdata/speed.RData"))
speed <- speed %>% 
  group_by(uuid, condition) %>% 
  mutate(total_time = max(end) - min(start), n=n()) %>% 
  mutate(bpm = n / (total_time / 60)) %>% 
  select(-n, -total_time) %>% 
  ungroup
speed <- speed %>% 
  mutate(uuid = paste("speed", uuid, sep=".")) %>% 
  mutate(condition = paste("speed", condition, sep="."))
save(speed, file = "data/rdata/speed2.RData")

# half_speed -------------------------------------------------------------------
# This creates a version of speed with only every odd item
# to compare against towse16
speed <- get(load("data/rdata/speed2.RData"))
half_speed <- speed %>% 
  filter(index %% 2 == 1) %>% 
  mutate(index = (index - 1) / 2 + 1)
half_speed <- half_speed %>% 
  mutate(across(c(uuid, condition), \(x){str_replace(x, "speed", "half_speed")}))
half_speed <- half_speed %>% 
  mutate(bpm = bpm / 2)
save(half_speed, file = "data/rdata/half_speed.RData")
