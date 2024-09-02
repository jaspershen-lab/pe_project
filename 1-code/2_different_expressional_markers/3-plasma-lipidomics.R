library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')
library(tidyverse)

load(
  "3-data_analysis/1-data-preparation/4-plasma-lipidomics-data/plasma_lipidomics_data.RData"
)

dir.create(
  "3-data_analysis/2_different_expressional_markers/3-plasma-lipidomics-data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/2_different_expressional_markers/3-plasma-lipidomics-data")

plasma_lipidomics_data@sample_info %>%
  dplyr::count(Trimester)

plasma_lipidomics_data@sample_info %>%
  ggplot(aes(GA_day, subject_id)) +
  geom_point(aes(color = class))

###Trimester 1
t1_data <-
  plasma_lipidomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester == "T1")

t1_data@sample_info %>%
  count(subject_id, Trimester)

t1_data@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()

###pca
pca_object <-
  t1_data %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t1_data, pca_object, color_by = "class")
plot
ggsave(plot,
       filename = "t1_pca_score_plot.pdf",
       width = 7,
       height = 6)


###p_value
control_sample_id <-
  t1_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "TC") %>%
  pull(sample_id)

pe_sample_id <-
  t1_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "PE") %>%
  pull(sample_id)

t1_data <-
  t1_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(10) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot <-
  t1_data %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    labs_x = "log2(Fold change, PE/TC)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_cutoff = 0.05,
    add_text = TRUE,
    text_from = "variable_id"
  )

plot

ggsave(plot,
       filename = "t1_volcano_plot.pdf",
       width = 7,
       height = 6)





###Trimester 2 early
t2_data_early <-
  plasma_lipidomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester == "T2" & GA_day < 150)

t2_data_early@sample_info %>%
  count(subject_id, Trimester)

t2_data_early@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()


###pca
pca_object <-
  t2_data_early %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_early, pca_object, color_by = "class")
plot
ggsave(plot,
       filename = "t2_early_pca_score_plot.pdf",
       width = 7,
       height = 6)


###p_value
control_sample_id <-
  t2_data_early %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "TC") %>%
  pull(sample_id)

pe_sample_id <-
  t2_data_early %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "PE") %>%
  pull(sample_id)

t2_data_early <-
  t2_data_early %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(10) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot <-
  t2_data_early %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    labs_x = "log2(Fold change, PE/TC)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_cutoff = 0.05,
    add_text = TRUE,
    text_from = "variable_id"
  )

plot

ggsave(plot,
       filename = "t2_early_volcano_plot.pdf",
       width = 7,
       height = 6)





###Trimester 2 late
t2_data_late <-
  plasma_lipidomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester == "T2" & GA_day > 150)

t2_data_late@sample_info %>%
  count(subject_id, Trimester)

t2_data_late@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()


###pca
pca_object <-
  t2_data_late %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_late, pca_object, color_by = "class")
plot
ggsave(plot,
       filename = "t2_late_pca_score_plot.pdf",
       width = 7,
       height = 6)

###p_value
control_sample_id <-
  t2_data_late %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "TC") %>%
  pull(sample_id)

pe_sample_id <-
  t2_data_late %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "PE") %>%
  pull(sample_id)

t2_data_late <-
  t2_data_late %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(10) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot <-
  t2_data_late %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    labs_x = "log2(Fold change, PE/TC)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_cutoff = 0.05,
    add_text = TRUE,
    text_from = "variable_id"
  )

plot

ggsave(plot,
       filename = "t2_late_volcano_plot.pdf",
       width = 7,
       height = 6)





###Trimester 3
t3_data <-
  plasma_lipidomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester == "T3")

t3_data@sample_info %>%
  count(subject_id, Trimester)

t3_data@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()


###pca
pca_object <-
  t3_data %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t3_data, pca_object, color_by = "class")
plot
ggsave(plot,
       filename = "t3_pca_score_plot.pdf",
       width = 7,
       height = 6)


###p_value
control_sample_id <-
  t3_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "TC") %>%
  pull(sample_id)

pe_sample_id <-
  t3_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "PE") %>%
  pull(sample_id)

t3_data <-
  t3_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(10) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot <-
  t3_data %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    labs_x = "log2(Fold change, PE/TC)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_cutoff = 0.05,
    add_text = TRUE,
    text_from = "variable_id"
  )

plot

ggsave(plot,
       filename = "t3_volcano_plot.pdf",
       width = 7,
       height = 6)
