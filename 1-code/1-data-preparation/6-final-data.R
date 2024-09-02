library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')
library(tidyverse)

####plasma metabolomics data
load(
  "3-data_analysis/1-data-preparation/2-plasma-metabolomics-data/metabolite/plasma_metabolomics_data.RData"
)

###urine metabolomics data
load(
  "3-data_analysis/1-data-preparation/3-urine-metabolomics-data/metabolite/urine_metabolomics_data.RData"
)

###plasma lipidomics data
load(
  "3-data_analysis/1-data-preparation/4-plasma-lipidomics-data/plasma_lipidomics_data.RData"
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
