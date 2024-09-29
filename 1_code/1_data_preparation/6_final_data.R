library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

####plasma metabolomics data
load(
  "3_data_analysis/1_data_preparation/2_plasma_metabolomics_data/metabolite/plasma_metabolomics_data.RData"
)

###urine metabolomics data
load(
  "3_data_analysis/1_data_preparation/3-urine-metabolomics-data/metabolite/urine_metabolomics_data.RData"
)

###plasma lipidomics data
load(
  "3_data_analysis/1_data_preparation/4-plasma-lipidomics-data/plasma_lipidomics_data.RData"
)

dim(plasma_metabolomics_data)
dim(plasma_lipidomics_data)
dim(urine_metabolomics_data)

colnames(plasma_metabolomics_data)
colnames(plasma_lipidomics_data)
colnames(urine_metabolomics_data)

length(intersect(
  colnames(plasma_metabolomics_data),
  colnames(plasma_lipidomics_data)
))
length(intersect(
  colnames(plasma_metabolomics_data),
  colnames(urine_metabolomics_data)
))

setdiff(colnames(plasma_metabolomics_data),
        colnames(plasma_lipidomics_data))
setdiff(colnames(plasma_lipidomics_data),
        colnames(plasma_metabolomics_data))
setdiff(colnames(urine_metabolomics_data),
        colnames(plasma_metabolomics_data))
