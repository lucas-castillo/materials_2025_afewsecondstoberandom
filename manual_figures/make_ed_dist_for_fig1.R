library(tidyverse)
# Load a distribution of empirical p values from a random sequence
load("output/ratio-measures/cooper12/processed/eps/cooper12.11__cooper12.a__1.RData")
temp_iid.distance %>% 
  ggplot(aes(ed)) + 
  geom_histogram(bins = 80, fill="#99c1f1") + 
  theme_void() + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))

ggsave("manual_figures/ed_dist.png", width = 12, height = 7)  
## Posterior distribution
set.seed(23)
x <- c("NO", 1:7, "H")
y <- dpois(1:9, 5)
tibble(x, y) %>% 
  mutate(x = factor(x, levels=x)) %>% 
  ggplot(aes(x,y)) +
  geom_col(fill="#1ca56f", color="black") + 
  theme_void() + 
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0))
ggsave("manual_figures/forest_allocation.png", width = 12, height = 7)  
