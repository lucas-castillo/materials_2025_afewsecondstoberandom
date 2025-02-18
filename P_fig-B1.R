library(tidyverse)
library(glue)
library(latex2exp)
library(patchwork)
# TR ----------------------------------------------------------------------
plot_tr <- function(df, title="", color_idx=2){
  identity <- TeX("$f(x) = 1x$")
  half_identity <- TeX("$f(x) = \\frac{1}{2}x$")
  
  p <- df %>% 
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
      fill=RColorBrewer::brewer.pal(3, "Set1")[color_idx]
    ) +
    geom_function(fun=\(x){x}) + 
    geom_function(fun=\(x){1/2*x}, linetype="dashed") + 
    scale_x_continuous(breaks=0:10, limits=c(1,7), expand=c(.035,0)) + 
    scale_y_continuous(breaks=0:10, limits=c(1,7), expand=c(.035,0)) + 
    annotate("text", x = 4, y = 4.5, label=identity, angle=55, size=5) +
    annotate("text", x = 5, y = 2.9, label=half_identity, angle=70/2, size=5) + 
    xlab("Thinning Rate (Fast)") +
    ylab("Thinning Rate (Slow)") + 
    theme_minimal(20) + 
    labs(title = title)
  p
}

plots <- list(
  plot_tr(readRDS("data/plot_data/speedtimes_speed2.rds"), title = "A. Our Experiment",2),
  plot_tr(readRDS("data/plot_data/speedtimes_towse16.rds"), title = "B. Towse et al. 2016",1),
  plot_tr(readRDS("data/plot_data/speedtimes_half_speed.rds"), title = "C. Our Experiment, halved",3)
)

p <- wrap_plots(plots)


ggsave(plot = p, "output/plots/figB1.jpg", dpi=900, width=13.9, height=5.87)
