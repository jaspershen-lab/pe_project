no_souce()

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

load(
  "3_data_analysis/1_data_preparation/2_plasma_metabolomics_data/metabolite/plasma_metabolomics_data.RData"
)

load(
  "3_data_analysis/1_data_preparation/3-urine-metabolomics-data/metabolite/urine_metabolomics_data.RData"
)

load(
  "3_data_analysis/1_data_preparation/4-plasma-lipidomics-data/plasma_lipidomics_data.RData"
)

dir.create("3_data_analysis/3_summary/",
           recursive = TRUE,
           showWarnings = FALSE)

setwd("3_data_analysis/3_summary/")


plasma_metabolomics_data@sample_info %>%
  ggplot(aes(GA_day, subject_id)) +
  geom_point(aes(color = Trimester_New))


plasma_lipidomics_data@sample_info %>%
  ggplot(aes(GA_day, subject_id)) +
  geom_point(aes(color = Trimester_New))

urine_metabolomics_data@sample_info %>%
  ggplot(aes(GA_day, subject_id)) +
  geom_point(aes(color = Trimester_New))


library(ggstatsplot)

sample_info <-
  plasma_lipidomics_data@sample_info %>%
  dplyr::distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::mutate(class = factor(class, levels = c("TC", "PE")))

sample_info %>%
  ggbetweenstats(class, Mom.Age) +
  theme_base +
  labs(x = "")

sample_info %>%
  ggbetweenstats(class, CR_length) +
  theme_base +
  labs(x = "")

sample_info %>%
  ggbetweenstats(class, Baby_weight) +
  theme_base +
  labs(x = "")

sample_info %>%
  ggbetweenstats(class, BMI_prePregnant) +
  theme_base +
  labs(x = "")

sample_info %>%
  group_by(class, Race) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(class) %>%
  dplyr::mutate(perc = n * 100 / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = class, y = perc)) +
  geom_bar(aes(fill = Race), position = "stack", stat = "identity") +
  theme_base +
  labs(x = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

sample_info %>%
  group_by(class, Ethnicity) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(class) %>%
  dplyr::mutate(perc = n * 100 / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = class, y = perc)) +
  geom_bar(aes(fill = Ethnicity), position = "stack", stat = "identity") +
  theme_base +
  labs(x = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

sample_info %>%
  group_by(class, Baby_sex) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(class) %>%
  dplyr::mutate(perc = n * 100 / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = class, y = perc)) +
  geom_bar(aes(fill = Baby_sex), position = "stack", stat = "identity") +
  theme_base +
  labs(x = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))





