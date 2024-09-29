###we don't have the urine lipidomics data

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidyverse)

load("3_data_analysis/1_data_preparation/1-phenotype-data/phenotype_data.RData")
