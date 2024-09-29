library(tidyverse)
library(ggplot2)

disease_color <-
  c("TC" = "blue",
    "PE" = "orange")

marker_color <-
  c("UP" = "#EE0000FF",
    "DOWN" = "#3B4992FF",
    "NO" = "#808180FF")


theme_base <-
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))