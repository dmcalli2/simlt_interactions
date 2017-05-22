#02 Arrange into ragged arrays
source("01_simulate_outcome_data.r")

library(glmmBUGS)

# Rename and convert to numeric
mydf <- filter(sim_data, drug_group_allocation == drugs_select) %>%
  arrange(desc(my_condition), desc(my_drug), studyno, strata, alloc)
mydf$study <- as.numeric(factor(mydf$studyno, levels = unique(mydf$studyno), labels = 1:length(unique(mydf$studyno))))
mydf$studyno <- NULL
mydf$events <- mydf$x
mydf$positive <- mydf$strata
mydf$indication <-as.numeric(factor(mydf$my_drug, levels = unique(mydf$my_drug), labels = 1:length(unique(mydf$my_drug))))
mydf <- as.data.frame (mydf)

# Create categorical and continuous data analysis ragged arrays
myrag_cat <- winBugsRaggedArray(mydf, effects = c("indication", "study"), observations = "events", 
                            covariates  = list(observations = c("positive", "alloc", "n")))

myrag_con <- winBugsRaggedArray(mydf, effects = c("indication", "study"), observations = "y", 
                                covariates  = list(observations = c("positive", "alloc", "se_prec")))
