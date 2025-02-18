library(tidyverse)
library(glue)
library(latex2exp)
library(patchwork)

# SPI ---------------------------------------------------------------------
voice <- 4.25
no_voice <- 2.30
pooled_df <- readRDS("data/plot_data/spi_pooled_df.rds")
names <- sort(unique(pooled_df$condition))
fancy_names <- c(
  "(1) Ca24.1N",
  "(1) Ca24.1U",
  "(2) Ca24.2",
  "(3) Co12", 
  "(4) Gr22",
  "(5) Gu23.ER",
  "(6) Gu23.IR",
  "(7) Gu23.MT",
  "(11) CE.F",
  "(11) CE.S",
  "(8) To16.1F", 
  "(8) To16.1S", 
  "(9) To16.2F", 
  "(9) To16.2S", 
  "(10) To16.3F", 
  "(10) To16.3S" 
)
fancy_names <- c(
  "(1) Castillo24 Exp. 1 Gauss.",
  "(1) Castillo24 Exp. 1 Unif.",
  "(2) Castillo24 Exp. 2",
  "(3) Cooper12", 
  "(4) Groot22",
  "(5) Guseva23 ER",
  "(6) Guseva23 IR",
  "(7) Guseva23 MT",
  "(11) Curr. Study Fast",
  "(11) Curr. Study Slow",
  "(8) Towse16 Exp. 1 Fast", 
  "(8) Towse16 Exp. 1 Slow", 
  "(9) Towse16 Exp. 2 Fast", 
  "(9) Towse16 Exp. 2 Slow", 
  "(10) Towse16 Exp. 3 Fast",
  "(10) Towse16 Exp. 3 Slow"
)
p1 <- pooled_df %>%
  mutate(condition = sapply(condition, \(x){fancy_names[which(names == x)]})) %>% 
  ungroup %>% 
  mutate(bpm = ifelse(BPM <= 30, "<= 30", ifelse(BPM > 50, "> 50", "30-50"))) %>% 
  mutate(bpm = factor(bpm, levels=c("<= 30", "30-50", "> 50"))) %>% 
  mutate(lb = Est - qnorm(0.975) * SE) %>%
  mutate(ub = Est + qnorm(0.975) * SE) %>% 
  ggplot(aes(forcats::fct_reorder(condition, -Est), Est)) + 
  geom_segment(aes(y=lb, yend=ub)) + 
  geom_point(aes(fill=domain, shape=output, size=bpm), stroke=.75) + 
  geom_hline(yintercept = voice, linetype="dashed") + 
  geom_hline(yintercept = no_voice, linetype="dashed") + 
  scale_y_continuous(
    breaks=c(0:10* 2),
    minor_breaks = NULL, 
    limits = c(0,10)
  ) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  scale_shape_manual(values = 21:25) + 
  scale_size_manual(values=2:4 / 3 * 4) + 
  ylab("Estimated seconds per item") +
  xlab("Dataset") +
  coord_flip() + 
  
  annotate("text", size=4.5, x = 4.5, y = no_voice - .35, angle=90, label=glue("Keypress or Click = {format(no_voice, nsmall=2)} spi")) +
  annotate("text", size=4.5, x = 13, y = voice - .35, angle=90, label=glue("Voice = {format(voice, nsmall=2)} spi")) +
  
  theme_minimal(20) + 
  theme(
    legend.box.background = element_rect(fill="#ffffffbb"),
    legend.position = c(.8, .7),
    legend.spacing.y = unit(-.5, "cm"),
    legend.title = element_text(margin = margin(b=5))
  ) + 
  guides(
    size=guide_legend(title="IPM", override.aes = list(fill="white", shape=21), byrow=T, order = 3),
    shape=guide_legend(title="Output Type", override.aes = list(size=4, fill="white"), byrow=T, order=2),
    fill=guide_legend(title="Domain", override.aes = list(shape=21, size=4), byrow=T, order = 1),
  )

p1
# TR ----------------------------------------------------------------------
df <- readRDS("data/plot_data/speedtimes_speed2.rds")
identity <- TeX("$f(x) = 1x$")
half_identity <- TeX("$f(x) = \\frac{1}{2}x$")

p2 <- df %>% 
  pivot_wider(
    names_from = condition, 
    values_from = rp, 
    id_cols = c(uuid, fastFirst)
  ) %>% 
  drop_na() %>% 
  ggplot(aes(Fast, Slow)) + 
  geom_point(
    size=4, 
    alpha=.6, 
    shape=21, 
    stroke=1,
    color="black", 
    fill=RColorBrewer::brewer.pal(3, "Set1")[2]
  ) +
  geom_function(fun=\(x){x}) + 
  geom_function(fun=\(x){1/2*x}, linetype="dashed") + 
  scale_x_continuous(breaks=0:10, limits=c(1,7), expand=c(.035,0)) + 
  scale_y_continuous(breaks=0:10, limits=c(1,7), expand=c(.035,0)) + 
  annotate("text", x = 4, y = 4.5, label=identity, angle=55, size=5) +
  annotate("text", x = 5, y = 2.9, label=half_identity, angle=70/2, size=5) + 
  xlab("Thinning Rate (Fast)") +
  ylab("Thinning Rate (Slow)") + 
  theme_minimal(20)

p1 + p2 + plot_annotation(tag_levels = "A")
ggsave("output/plots/fig2.jpg", dpi=600, width=11.3*1.1, height=6.13*1.1)

