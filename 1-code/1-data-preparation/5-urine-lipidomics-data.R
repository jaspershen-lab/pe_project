library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')
library(tidyverse)

load("3-data_analysis/1-data-preparation/1-phenotype-data/phenotype_data.RData")
