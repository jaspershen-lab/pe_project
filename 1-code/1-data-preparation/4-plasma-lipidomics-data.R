library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')
library(tidyverse)

load("3-data_analysis/1-data-preparation/1-phenotype-data/phenotype_data.RData")

data <-
  readr::read_csv("2-data/Lipidomics_Preeclampsia_data.csv")

dir.create(
  "3-data_analysis/1-data-preparation/4-plasma-lipidomics-data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/1-data-preparation/4-plasma-lipidomics-data/")

phenotype_data$URINE_ID
phenotype_data$Plasma_ID

phenotype_data$sample_id_urine
phenotype_data$sample_id_plasma

phenotype_data$Sample_ID_urine
phenotype_data$Sample_ID_plasma

phenotype_data$Patient_ID

match(colnames(data), paste0("X", phenotype_data$Plasma_ID))

colnames(data)[1] <-
  "variable_id"

expression_data <-
  data %>%
  tibble::column_to_rownames(var = "variable_id")

sample_info <-
  data.frame(sample_id = colnames(expression_data))

variable_info <-
  data.frame(variable_id = rownames(expression_data))

phenotype_data <-
  phenotype_data %>%
  dplyr::select(Plasma_ID, URINE_ID, Measurement_order:sample_id_universal)

phenotype_data$Plasma_ID <-
  paste("X", phenotype_data$Plasma_ID, sep = "")

phenotype_data$URINE_ID <-
  paste("X", phenotype_data$URINE_ID, sep = "")

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = c("sample_id" = "Plasma_ID"))

sample_info <-
  sample_info %>%
  dplyr::rename(sample_id_plasma_lipidomics = sample_id,
                sample_id = sample_id_universal) %>%
  dplyr::select(sample_id, dplyr::everything())

colnames(expression_data) <-
  sample_info$sample_id

library(tidymass)

plasma_lipidomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

save(plasma_lipidomics_data, file = "plasma_lipidomics_data.RData")
