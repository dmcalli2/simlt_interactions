# 00 Read in data
library(tidyverse)

## Read in diabetes trials
load("data/metadata_for_simulation.Rdata")


diabetes_final <- diabetes_final %>% 
  filter(atc_5 != "A10BF")

