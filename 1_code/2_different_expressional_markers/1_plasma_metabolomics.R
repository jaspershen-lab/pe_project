no_souce()

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)
library(tidymass)

load(
  "3_data_analysis/1_data_preparation/2_plasma_metabolomics_data/metabolite/plasma_metabolomics_data.RData"
)

dir.create(
  "3_data_analysis/2_different_expressional_markers/1_plasma_metabolomics_data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/2_different_expressional_markers/1_plasma_metabolomics_data")

plasma_metabolomics_data@sample_info$BMI_prePregnant[is.na(plasma_metabolomics_data@sample_info$BMI_prePregnant)] <-
  mean(plasma_metabolomics_data@sample_info$BMI_prePregnant[!is.na(plasma_metabolomics_data@sample_info$BMI_prePregnant)])

plasma_metabolomics_data@sample_info %>%
  dplyr::count(Trimester_New)

plasma_metabolomics_data@sample_info %>%
  ggplot(aes(GA_day, subject_id)) +
  geom_point(aes(color = class))

plasma_metabolomics_data@sample_info %>%
  count(Trimester_New)

plasma_metabolomics_data <-
  plasma_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::arrange(desc(class))

###Trimester 1
t1_data <-
  plasma_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester_New == "T1")

t1_data@sample_info %>%
  count(subject_id, Trimester_New)

t1_data@sample_info %>%
  ggplot(aes(x = 1, y = GA_day)) +
  geom_boxplot() +
  geom_jitter(aes(color = class), width = 0.2)

t1_data@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()

###p_value and fc
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

x <-
  t1_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

y <-
  t1_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(2) %>%
  adjust_confounder(confounder_name_list = c("BMI_prePregnant", "Race")) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot(-log(x@variable_info$p_value, 10),-log(y@variable_info$p_value, 10))
plot(
  -log(x@variable_info$p_value_adjust, 10),-log(y@variable_info$p_value_adjust, 10)
)
sum(x@variable_info$p_value_adjust < 0.05)
sum(y@variable_info$p_value_adjust < 0.05)

which(x@variable_info$p_value_adjust < 0.05)
which(y@variable_info$p_value_adjust < 0.05)

i = 440
plot(as.numeric(x@expression_data[i, ]),
     as.numeric(y@expression_data[i, ]))

wilcox.test(as.numeric(x@expression_data[i, control_sample_id]),
            as.numeric(x@expression_data[i, pe_sample_id]))

wilcox.test(as.numeric(y@expression_data[i, control_sample_id]),
            as.numeric(y@expression_data[i, pe_sample_id]))

t1_data <-
  t1_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(2) %>%
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
    text_from = "Compound_name"
  )

plot

ggsave(plot,
       filename = "t1_volcano_plot.pdf",
       width = 7,
       height = 6)

###pca
####all features
pca_object <-
  t1_data %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t1_data, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t1_pca_score_plot.pdf",
       width = 7,
       height = 6)

####all markers
pca_object <-
  t1_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::slice_head(n = 10) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t1_data, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t1_pca_score_plot_markers.pdf",
       width = 7,
       height = 6)

###heatmap
library(ComplexHeatmap)
###all features
# Create a block annotation
ha <- HeatmapAnnotation(
  Group = t1_data@sample_info$class,
  col = list(Group = c(
    "TC" = unname(disease_color["TC"]),
    "PE" = unname(disease_color["PE"])
  )),
  # Colors for the groups
  annotation_name_side = "left"  # Position of the group labels
)

plot <-
  t1_data %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t1_data@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t1_heatmap.pdf",
       width = 7,
       height = 6)


####all markers
plot <-
  t1_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::slice_head(n = 10) %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t1_data@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t1_heatmap_markers.pdf",
       width = 7,
       height = 6)

###Trimester 2 early
t2_data_early <-
  plasma_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester_New == "T2" & GA_day < 150)

t2_data_early@sample_info %>%
  count(subject_id, Trimester_New)

sort(t2_data_early@sample_info$GA_day)

t2_data_early@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()

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
  log(2) %>%
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
    text_from = "Compound_name"
  )

plot

ggsave(plot,
       filename = "t2_early_volcano_plot.pdf",
       width = 7,
       height = 6)

###pca
####all features
pca_object <-
  t2_data_early %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_early, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t2_early_pca_score_plot.pdf",
       width = 7,
       height = 6)

####all markers
pca_object <-
  t2_data_early %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_early, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t2_early_pca_score_plot_markers.pdf",
       width = 7,
       height = 6)

###heatmap
library(ComplexHeatmap)
###all features
# Create a block annotation
ha <- HeatmapAnnotation(
  Group = t2_data_early@sample_info$class,
  col = list(Group = c(
    "TC" = unname(disease_color["TC"]),
    "PE" = unname(disease_color["PE"])
  )),
  # Colors for the groups
  annotation_name_side = "left"  # Position of the group labels
)

plot <-
  t2_data_early %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t2_data_early@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t2_early_heatmap.pdf",
       width = 7,
       height = 6)


####all markers
plot <-
  t2_data_early %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t2_data_early@sample_info$class, levels = c("TC", "PE")),
    row_km = 2
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t2_early_heatmap_markers.pdf",
       width = 7,
       height = 6)



###Trimester 2 late
t2_data_late <-
  plasma_metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Trimester_New == "T2" & GA_day > 150)

t2_data_late@sample_info %>%
  count(subject_id, Trimester_New)

sort(t2_data_late@sample_info$GA_day)


t2_data_late@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()


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
  log(2) %>%
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
    text_from = "Compound_name"
  )

plot

ggsave(plot,
       filename = "t2_late_volcano_plot.pdf",
       width = 7,
       height = 6)

###pca
####all features
pca_object <-
  t2_data_late %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_late, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t2_late_pca_score_plot.pdf",
       width = 7,
       height = 6)

####all markers
pca_object <-
  t2_data_late %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(t2_data_late, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "t2_late_pca_score_plot_markers.pdf",
       width = 7,
       height = 6)

###heatmap
library(ComplexHeatmap)
###all features
# Create a block annotation
ha <- HeatmapAnnotation(
  Group = t2_data_late@sample_info$class,
  col = list(Group = c(
    "TC" = unname(disease_color["TC"]),
    "PE" = unname(disease_color["PE"])
  )),
  # Colors for the groups
  annotation_name_side = "left"  # Position of the group labels
)

plot <-
  t2_data_late %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t2_data_late@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t2_late_heatmap.pdf",
       width = 7,
       height = 6)


####all markers
plot <-
  t2_data_late %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(t2_data_late@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "t2_late_heatmap_markers.pdf",
       width = 7,
       height = 6)






###combine three trimesters
combined_data <-
  plasma_metabolomics_data %>%
  dplyr::filter(Trimester_New %in% c("T1", "T2")) %>%
  massdataset::summarise_samples(group_by = "subject_id", what = "median_intensity")

combined_data@sample_info %>%
  count(subject_id, Trimester_New)

sort(combined_data@sample_info$GA_day)


combined_data@expression_data[4, ] %>%
  as.numeric() %>%
  `+`(1) %>%
  log(2) %>%
  density() %>%
  plot()


###p_value
control_sample_id <-
  combined_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "TC") %>%
  pull(sample_id)

pe_sample_id <-
  combined_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "PE") %>%
  pull(sample_id)

combined_data <-
  combined_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    mean_median = "median"
  ) %>%
  activate_mass_dataset(what = "expression_data") %>%
  `+`(1) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = pe_sample_id,
    method = "wilcox.test",
    p_adjust_methods = "fdr"
  )

plot <-
  combined_data %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    labs_x = "log2(Fold change, PE/TC)",
    labs_y = "-log(p-adjust, 10)",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_cutoff = 0.05,
    add_text = TRUE,
    text_from = "Compound_name"
  )

plot

ggsave(plot,
       filename = "combined_volcano_plot.pdf",
       width = 7,
       height = 6)

###pca
####all features
pca_object <-
  combined_data %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(combined_data, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "combined_pca_score_plot.pdf",
       width = 7,
       height = 6)

####all markers
pca_object <-
  combined_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale() %>%
  run_pca()

plot <-
  pca_score_plot(combined_data, pca_object, color_by = "class") +
  scale_color_manual(values = disease_color) +
  scale_fill_manual(values = disease_color)

plot

ggsave(plot,
       filename = "combined_pca_score_plot_markers.pdf",
       width = 7,
       height = 6)


###heatmap
library(ComplexHeatmap)
###all features
# Create a block annotation
ha <- HeatmapAnnotation(
  Group = combined_data@sample_info$class,
  col = list(Group = c(
    "TC" = unname(disease_color["TC"]),
    "PE" = unname(disease_color["PE"])
  )),
  # Colors for the groups
  annotation_name_side = "left"  # Position of the group labels
)

plot <-
  combined_data %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(combined_data@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "combined_heatmap.pdf",
       width = 7,
       height = 6)


####all markers
plot <-
  combined_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::arrange(p_value_adjust) %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale_data() %>%
  extract_expression_data() %>%
  Heatmap(
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    name = "Z-score",
    top_annotation = ha,
    border = TRUE,
    column_split = factor(combined_data@sample_info$class, levels = c("TC", "PE"))
  )

plot <- ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "combined_heatmap_markers.pdf",
       width = 7,
       height = 6)


######markers in all trimesters
temp_data1 <-
  t1_data@variable_info %>%
  dplyr::select(variable_id, Compound_name, fc, p_value_adjust) %>%
  dplyr::mutate(trimester = "T1")

temp_data2 <-
  t2_data_early@variable_info %>%
  dplyr::select(variable_id, fc, p_value_adjust) %>%
  dplyr::mutate(trimester = "T2_early")

temp_data3 <-
  t2_data_late@variable_info %>%
  dplyr::select(variable_id, fc, p_value_adjust) %>%
  dplyr::mutate(trimester = "T2_late")

temp_data4 <-
  combined_data@variable_info %>%
  dplyr::select(variable_id, fc, p_value_adjust) %>%
  dplyr::mutate(trimester = "Combined")

temp_data <-
  temp_data1 %>%
  full_join(temp_data2,
            by = c("variable_id", "fc", "p_value_adjust", "trimester")) %>%
  full_join(temp_data3,
            by = c("variable_id", "fc", "p_value_adjust", "trimester")) %>%
  full_join(temp_data4,
            by = c("variable_id", "fc", "p_value_adjust", "trimester")) %>%
  dplyr::mutate(
    up_down = case_when(
      fc > 1 & p_value_adjust < 0.05 ~ "UP",
      fc < 1 &
        p_value_adjust < 0.05 ~ "DOWN",
      TRUE ~ "NO"
    )
  ) %>%
  dplyr::mutate(up_down = factor(up_down, levels = c("NO", "UP", "DOWN"))) %>%
  dplyr::arrange(up_down) %>%
  dplyr::mutate(trimester = factor(trimester, levels = c("T1", "T2_early", "T2_late", "Combined")))

library(ggforce)

####use ggraph network to show the markers in all trimesters
node_data <-
  temp_data %>%
  dplyr::mutate(variable_id = paste(variable_id, trimester, sep = "_"))

edge_data <-
  rbind(
    data.frame(
      from = paste(temp_data$variable_id, "T1", sep = "_"),
      to = paste(temp_data$variable_id, "T2_early", sep = "_")
    ),
    data.frame(
      from = paste(temp_data$variable_id, "T2_early", sep = "_"),
      to = paste(temp_data$variable_id, "T2_late", sep = "_")
    ),
    data.frame(
      from = paste(temp_data$variable_id, "T2_late", sep = "_"),
      to = paste(temp_data$variable_id, "Combined", sep = "_")
    )
  ) %>%
  dplyr::left_join(node_data[, c("variable_id", "fc")], by = c("from" = "variable_id")) %>%
  dplyr::rename(fc_from = fc) %>%
  dplyr::left_join(node_data[, c("variable_id", "fc")], by = c("to" = "variable_id")) %>%
  dplyr::rename(fc_to = fc)

edge_data <-
  edge_data %>%
  dplyr::mutate(
    edge_class = case_when(
      fc_from > 1 & fc_to > 1 ~ "same",
      fc_from < 1 & fc_to < 1 ~ "same",
      fc_from > 1 & fc_to < 1 ~ "different",
      fc_from < 1 & fc_to > 1 ~ "different",
      TRUE ~ "different"
    )
  )

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE)

g <- total_graph
igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
coords <-
  ggraph::create_layout(g, layout = "bipartite")

coords$x[coords$trimester == "T1"] = 1
coords$x[coords$trimester == "T2_early"] = 2
coords$x[coords$trimester == "T2_late"] = 3
coords$x[coords$trimester == "Combined"] = 4

coords$y <-
  log(coords$fc, 2)

library(ggraph)
library(tidygraph)
library(igraph)
my_graph <-
  ggraph::create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph::ggraph(my_graph, layout = 'bipartite') +
  ggraph::geom_edge_diagonal(
    strength = 0.5,
    aes(color = edge_class),
    edge_width = 0.5,
    alpha = 0.2,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    aes(fill = up_down),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = marker_color) +
  theme_base +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "", y = "log(FC, 2)") +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("T1", "T2_early", "T2_late", "Combined")
  )

plot
ggsave(plot,
       filename = "features_in_all_trimesters.pdf",
       width = 6,
       height = 8)


###remain the features significant in at least one trimester
marker_id <-
  temp_data %>%
  dplyr::filter(up_down != "NO") %>%
  pull(variable_id) %>%
  unique()

node_data <-
  temp_data %>%
  dplyr::filter(variable_id %in% marker_id) %>%
  dplyr::mutate(variable_id = paste(variable_id, trimester, sep = "_"))

edge_data <-
  rbind(
    data.frame(
      from = paste(temp_data$variable_id, "T1", sep = "_"),
      to = paste(temp_data$variable_id, "T2_early", sep = "_")
    ),
    data.frame(
      from = paste(temp_data$variable_id, "T2_early", sep = "_"),
      to = paste(temp_data$variable_id, "T2_late", sep = "_")
    ),
    data.frame(
      from = paste(temp_data$variable_id, "T2_late", sep = "_"),
      to = paste(temp_data$variable_id, "Combined", sep = "_")
    )
  ) %>%
  dplyr::left_join(node_data[, c("variable_id", "fc")], by = c("from" = "variable_id")) %>%
  dplyr::rename(fc_from = fc) %>%
  dplyr::left_join(node_data[, c("variable_id", "fc")], by = c("to" = "variable_id")) %>%
  dplyr::rename(fc_to = fc) %>%
  dplyr::filter(from %in% node_data$variable_id &
                  to %in% node_data$variable_id)

edge_data <-
  edge_data %>%
  dplyr::mutate(
    edge_class = case_when(
      fc_from > 1 & fc_to > 1 ~ "same",
      fc_from < 1 & fc_to < 1 ~ "same",
      fc_from > 1 & fc_to < 1 ~ "different",
      fc_from < 1 & fc_to > 1 ~ "different",
      TRUE ~ "different"
    )
  )

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE)


g <- total_graph
igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
coords <-
  ggraph::create_layout(g, layout = "bipartite")

coords$x[coords$trimester == "T1"] = 1
coords$x[coords$trimester == "T2_early"] = 2
coords$x[coords$trimester == "T2_late"] = 3
coords$x[coords$trimester == "Combined"] = 4

coords$y <-
  log(coords$fc, 2)

my_graph <-
  ggraph::create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  ) %>%
  dplyr::filter()

plot <-
  ggraph::ggraph(my_graph, layout = 'bipartite') +
  ggraph::geom_edge_diagonal(
    strength = 0.5,
    aes(color = edge_class),
    edge_width = 0.5,
    alpha = 0.2,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    aes(fill = up_down, size = -log(p_value_adjust, 10)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_fill_manual(values = marker_color) +
  theme_base +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "", y = "log(FC, 2)") +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("T1", "T2_early", "T2_late", "Combined")
  )
plot
ggsave(plot,
       filename = "markers_in_all_trimesters.pdf",
       width = 7,
       height = 8)
