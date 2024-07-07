library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')
library(tidyverse)

load("2-data/PE_urine_input_final.RData")
load("3-data_analysis/1-data-preparation/1-phenotype-data/phenotype_data.RData")


dir.create(
  "3-data_analysis/1-data-preparation/3-urine-metabolomics-data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3-data_analysis/1-data-preparation/3-urine-metabolomics-data/")

###feature table
expression_data <-
  urine_comp %>% 
  as.data.frame()

sample_info <-
  data.frame(sample_id = colnames(expression_data))

variable_info <-
  data.frame(variable_id = rownames(expression_data))

PE_annotation_urine <-
  PE_annotation_urine %>%
  dplyr::rename(variable_id = Compound)

variable_info <-
  variable_info %>%
  dplyr::left_join(PE_annotation_urine, by = "variable_id")

variable_info <-
  variable_info %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "\\/", "")) %>%
  dplyr::rename(Compound_name = Metabolite)

rownames(expression_data) <-
  variable_info$variable_id

phenotype_data <-
  phenotype_data %>%
  dplyr::select(sample_id_urine,
                sample_id_urine,
                Measurement_order:sample_id_universal)

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = c("sample_id" = "sample_id_urine"))

sample_info <-
  sample_info %>%
  dplyr::rename(sample_id_urine_metabolomics = sample_id,
                sample_id = sample_id_universal) %>%
  dplyr::select(sample_id, dplyr::everything())

sample_info$subject_id

colnames(expression_data) <-
  sample_info$sample_id

dir.create("feature")
library(tidymass)

urine_metabolomics_data <-
  massdataset::create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

save(urine_metabolomics_data, file = "feature/urine_metabolomics_data.RData")

###metabolite data
expression_data <-
  urine_annotated_data %>% 
  as.data.frame()

sample_info <-
  data.frame(sample_id = colnames(expression_data))

variable_info <-
  data.frame(variable_id = rownames(expression_data))

variable_info <-
  variable_info %>%
  dplyr::left_join(PE_annotation_urine, by = "variable_id")

variable_info <-
  variable_info %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "\\/", "")) %>%
  dplyr::rename(Compound_name = Metabolite)

rownames(expression_data) <-
  variable_info$variable_id

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = c("sample_id" = "sample_id_urine"))

sample_info <-
  sample_info %>%
  dplyr::rename(sample_id_urine_metabolomics = sample_id,
                sample_id = sample_id_universal) %>%
  dplyr::select(sample_id, dplyr::everything())

sample_info$subject_id

colnames(expression_data) <-
  sample_info$sample_id

dir.create("metabolite")
library(tidymass)

urine_metabolomics_data <-
  massdataset::create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

save(urine_metabolomics_data, file = "metabolite/urine_metabolomics_data.RData")
