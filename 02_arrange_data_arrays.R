#02 Arrange into ragged arrays
library(glmmBUGS)
library(tidyverse)

# Read in sim_study
sim_data <- readRDS(file = "scratch_data/simulated_data.Rds")

# Add study number
sim_data <- sim_data %>% 
  group_by(drug_group_allocation) %>% 
  mutate(studyno = 1:length(drug_group_allocation)) %>% 
  ungroup()

# Rename and convert to numeric
# drugs_select <- "Airways"

mydf <- filter(sim_data, drug_group_allocation == drugs_select) %>%
  arrange(desc(my_condition), desc(my_drug), studyno)
mydf$study <- as.numeric(factor(mydf$studyno, levels = unique(mydf$studyno),
                                labels = 1:length(unique(mydf$studyno))))
mydf$studyno <- NULL
mydf$drug <-as.numeric(factor(mydf$my_drug, levels = unique(mydf$my_drug),
                                    labels = 1:length(unique(mydf$my_drug))))
mydf$condition <-as.numeric(factor(mydf$my_condition, levels = unique(mydf$my_condition),
                                    labels = 1:length(unique(mydf$my_condition))))
mydf <- as.data.frame (mydf)

# Create ragged array
myrag <- winBugsRaggedArray(mydf, effects = c("drug", "condition"),
                            observations = c("con", "cat", "log"), 
                            covariates  = list(observations = c("study")))

