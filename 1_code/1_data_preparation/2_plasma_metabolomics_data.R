library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

load("2-data/PE_plasma_input_final.RData")
load("3_data_analysis/1_data_preparation/1-phenotype-data/phenotype_data.RData")


dir.create(
  "3_data_analysis/1_data_preparation/2_plasma_metabolomics_data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/1_data_preparation/2_plasma_metabolomics_data/")

###feature table
expression_data <-
  plasma_comp

sample_info <-
  data.frame(sample_id = colnames(expression_data))

variable_info <-
  data.frame(variable_id = rownames(expression_data))

PE_annotation_plasma <-
  PE_annotation_plasma %>%
  dplyr::rename(variable_id = Compound)

variable_info <-
  variable_info %>%
  dplyr::left_join(PE_annotation_plasma, by = "variable_id")

variable_info <-
  variable_info %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "\\/", "")) %>%
  dplyr::rename(Compound_name = Metabolite)

rownames(expression_data) <-
  variable_info$variable_id

phenotype_data <-
  phenotype_data %>%
  dplyr::select(sample_id_plasma,
                sample_id_urine,
                Measurement_order:sample_id_universal)

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = c("sample_id" = "sample_id_plasma"))

sample_info <-
  sample_info %>%
  dplyr::rename(sample_id_plasma_metabolomics = sample_id,
                sample_id = sample_id_universal) %>%
  dplyr::select(sample_id, dplyr::everything())

colnames(expression_data) <-
  sample_info$sample_id

dir.create("feature")
library(tidymass)

plasma_metabolomics_data <-
  massdataset::create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

save(plasma_metabolomics_data, file = "feature/plasma_metabolomics_data.RData")

###metabolite data
expression_data <-
  plasma_annotated_data

sample_info <-
  data.frame(sample_id = colnames(expression_data))

variable_info <-
  data.frame(variable_id = rownames(expression_data))

variable_info <-
  variable_info %>%
  dplyr::left_join(PE_annotation_plasma, by = "variable_id")

variable_info <-
  variable_info %>%
  dplyr::mutate(variable_id = stringr::str_replace(variable_id, "\\/", "")) %>%
  dplyr::rename(Compound_name = Metabolite)

rownames(expression_data) <-
  variable_info$variable_id

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = c("sample_id" = "sample_id_plasma"))

sample_info <-
  sample_info %>%
  dplyr::rename(sample_id_plasma_metabolomics = sample_id,
                sample_id = sample_id_universal) %>%
  dplyr::select(sample_id, dplyr::everything())

colnames(expression_data) <-
  sample_info$sample_id

dir.create("metabolite")
library(tidymass)

plasma_metabolomics_data <-
  massdataset::create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

save(plasma_metabolomics_data, file = "metabolite/plasma_metabolomics_data.RData")
