#01 Simulate outcome data)
load("Dataset_of_trials.Rdata")

# Packages
library(dplyr)
library(tidyr)

sim_data <- select (trials, drug_group_allocation, my_condition, my_drug, enrollment)
rm(trials)

## Categorical data constants
comorbidity_prevalence <- 0.2
baseline_risk <- 0.2
comorbidity_effect <- runif(nrow(sim_data), 1.2, 1.5)
tx_effect <- runif(nrow(sim_data), 0.65,0.75)
interaction_effect <- 1.0

## COntinuous data constants
baseline_mean <- 70
mean_tx_effect <- rnorm(nrow(sim_data), 10, 1)
mean_comorbidity_effect <- rnorm(nrow(sim_data), 5, 0.5)
mean_interaction_effect <- 0 #2.5
baseline_sd <- rnorm(nrow(sim_data), 15, 1)


sim_data$pl_n_0 <- sim_data$enrollment * 0.5 * (1-comorbidity_prevalence)
sim_data$pl_n_1 <- sim_data$enrollment * 0.5 * comorbidity_prevalence
sim_data$tx_n_0 <- sim_data$pl_n_0
sim_data$tx_n_1 <- sim_data$pl_n_1

## Events in placebo
sim_data$pl_x_0 <- sim_data$pl_n_0 * baseline_risk
sim_data$pl_x_1 <- sim_data$pl_n_1 * baseline_risk * comorbidity_effect

## Events in treatment
sim_data$tx_x_0 <- sim_data$tx_n_0 * baseline_risk * tx_effect
sim_data$tx_x_1 <- sim_data$tx_n_1 * baseline_risk * comorbidity_effect * tx_effect *interaction_effect

# Mean outcome in placebo
sim_data$pl_y_0 <- baseline_mean
sim_data$pl_y_1 <- baseline_mean + mean_comorbidity_effect

# Mean outcome in tx
sim_data$tx_y_0 <- baseline_mean - mean_tx_effect
sim_data$tx_y_1 <- baseline_mean + mean_comorbidity_effect - mean_tx_effect + mean_interaction_effect

# Round all
vars <- c('pl_n_0', 'pl_n_1', 'tx_n_0', 'tx_n_1', 'pl_x_0', 'pl_x_1', 'tx_x_0', 'tx_x_1')
sim_data [ , vars] <- lapply(sim_data[,vars], round, digits = 0)

# Assign study_number
sim_data$studyno <- seq_along(sim_data$drug_group_allocation)

## Move into long format needed for analysis
sim_data <- sim_data %>%
  select(-enrollment) %>%
  gather(varname, value, -drug_group_allocation, -my_drug, - my_condition, -studyno) %>%
  separate (varname, into = c("alloc", "event_n", "strata"), sep = "_") 
sim_data <- spread(sim_data, event_n, value)

## Assign standard error for each
sim_data$se_prec <- 1/(baseline_sd/ sim_data$n^0.5)^2

# Convert alloc and strata into numeric
sim_data$alloc <- ifelse(sim_data$alloc == "tx", 1, 0)
sim_data$strata <- as.numeric(sim_data$strata)


