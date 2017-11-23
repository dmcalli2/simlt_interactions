# 00 Read in data
library(tidyverse)

## Read in diabetes trials
load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")


diabetes_final <- diabetes_final %>% 
  filter(atc_5 != "A10BF")
