library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

# Load the data

data1 <-
  readr::read_csv("2-data/meta_data_lipidomics.csv")

data1 <-
  data1 %>%
  dplyr::rename(
    sample_id = Sample_ID,
    sample_id_plasma = Plasma_ID,
    subject_id = Patient_ID,
    injection_order_plasma = Measurement_order
  )

data2 <-
  readr::read_csv("2-data/metadata_plasma_urine_final.csv")

data2 <-
  data2 %>%
  dplyr::rename(
    sample_id_urine = Urine_sampleID,
    sample_id_plasma = Plasma_sampleID,
    subject_id = SubjectID,
    class = Patient_class
  )

phenotype_data <-
  data2

dir.create(
  "3_data_analysis/1_data_preparation/1-phenotype-data",
  recursive = TRUE,
  showWarnings = FALSE
)

setwd("3_data_analysis/1_data_preparation/1-phenotype-data/")

phenotype_data$URINE_ID <-
  as.character(phenotype_data$URINE_ID)

phenotype_data$Plasma_ID <-
  as.character(phenotype_data$Plasma_ID)

phenotype_data$subject_id <-
  as.character(phenotype_data$subject_id)

phenotype_data$Trimester_sub

temp <-
  phenotype_data %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = sum(duplicated(Trimester_sub)))

phenotype_data$sample_id_universal <-
  paste(phenotype_data$subject_id,
        phenotype_data$Trimester_sub,
        sep = "_")

save(phenotype_data, file = "phenotype_data.RData")
