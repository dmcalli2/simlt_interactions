# 03_summarise_simulations
library(tidyverse)
library(stringr)
library(ggplot2)

scenarios <- c("scenario_atc5_0.05_trial_0.05_drug_0.05.rds",
               "scenario_atc5_0.1_trial_0.1_drug_0.1.rds",
               "scenario_atc5_0.1_trial_0.1_drug_0.25.rds",
               "scenario_atc5_0.1_trial_0.25_drug_0.1.rds", 
               "scenario_atc5_0.25_trial_0.1_drug_0.1.rds",
               "scenario_atc5_0.25_trial_0.25_drug_0.25.rds")

scenarios_names <- rbind(
  c(0.05, 0.05, 0.05),
  c(0.1, 0.1, 0.1),
  c(0.1, 0.1, 0.25),
  c(0.1, 0.25, 0.1), 
  c(0.25, 0.1, 0.1),
  c(0.25, 0.25, 0.25))
scenarios_names <- as.data.frame(scenarios_names)
names(scenarios_names) <- c("atc5", "trial", "drug")

scenario_res <- lapply(scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("simulate_interactions/effect_0_point2_priors0point01/", each_scenario))
  fxd <- lapply(each_scenario, function(x) x$fixed)
  fxd <- do.call(rbind, fxd)
  fxd <- as.data.frame(fxd)
  sapply(fxd, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
})
names(scenario_res) <- scenarios

scenario_res <- lapply(1:6, function(i) {
  scenario_res <- as.data.frame(scenario_res[[i]])
  scenario_res <- cbind(scenario_res, scenarios_names[i,])
  scenario_res$qs <- paste0("q", c(2.5, 50, 97.5))
  scenario_res
})

scenario_res <- bind_rows(scenario_res)


scenario_res2 <- scenario_res %>% 
  gather(key = "var", value = "value", -atc5, -drug, -trial, -qs) %>% 
  spread(key = qs, value = value)

plot1 <- ggplot(scenario_res2 %>% 
                  filter(var == "0.5quant") %>% 
                  mutate(trial = paste0("Trial ", trial),
                         drug = paste0("Drug ", drug)),
                aes(x = factor(atc5), y = q50, ymin = q2.5, ymax = q97.5)) +
  geom_errorbar() +
  geom_point() +
  facet_grid(trial ~ drug) +
  scale_x_discrete("Drug class") +
  scale_y_continuous("Effect estimate")
plot1
