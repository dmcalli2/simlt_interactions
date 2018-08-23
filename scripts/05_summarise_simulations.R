# 03_summarise_simulations
library(tidyverse)
library(stringr)
library(ggplot2)

# Identify saved results
sim1_scenarios <- list.files("unix_results/sim1/", recursive=TRUE, pattern = "sim1")
sim1_oneclass_scenarios <- list.files("unix_results/sim1/", recursive=TRUE, pattern = "oneclass")
sim1_scenarios <- setdiff(sim1_scenarios,sim1_oneclass_scenarios)

sim2_scenarios <- list.files("unix_results/sim2/", recursive=TRUE, pattern = "sim2")
sim2_oneclass_scenarios <- list.files("unix_results/sim2/", recursive=TRUE, pattern = "oneclass")
sim2_scenarios <- setdiff(sim2_scenarios,sim2_oneclass_scenarios)


ScenarioNames <- function (scenarios, sim) {
  if(sim==1){ scenarios_names <-  str_split_fixed( str_sub(scenarios, -37, -5), "_", n = 7)
  scenarios_names <- cbind(scenarios_names[,1],apply(scenarios_names[, c(3, 5, 7)], 2, as.double))
  scenarios_names <- as.data.frame(scenarios_names)
  names(scenarios_names) <- c("como_prev","atc5/moa", "trial", "drug")
  }
  else if(sim == 2) {scenarios_names <-  str_split_fixed( str_sub(scenarios, -46, -5), "_", n = 9)
  scenarios_names <- cbind(scenarios_names[,1],apply(scenarios_names[, c(3, 5, 7, 9)], 2, as.double))
  scenarios_names <- as.data.frame(scenarios_names)
  names(scenarios_names) <- c("como_prev","path","atc5/moa", "trial", "drug")
  }
  scenarios_names
}

sim1_scenarios_names <- s1ScenarioNames(sim1_scenarios,sim=1)
sim2_scenarios_names <- s1ScenarioNames(sim2_scenarios,sim=2)

# read and convert to data frame for each scenario
sim1_scenario_res <- lapply(sim1_scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("unix_results/sim1/", each_scenario))})

sim1_scenario_res <- map(sim1_scenario_res, function (scen){
  fxd <- map(scen, ~ .x$fixed) 
  fxd <- do.call(rbind, fxd)
  fxd <- as.tibble(fxd) %>% 
  mutate(iter = 1:1000)
  fxd
})
names(sim1_scenario_res) <- str_sub(sim1_scenarios, -37, -5)
sim1_scenario_res <- bind_rows(sim1_scenario_res, .id = "scenario")

# read and convert to data frame for each scenario
sim2_scenario_res <- lapply(sim2_scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("unix_results/sim2/", each_scenario))})

sim2_scenario_res <- map(sim2_scenario_res, function (scen){
  fxd <- map(scen, ~ .x$fixed) 
  fxd <- do.call(rbind, fxd)
  fxd <- as.tibble(fxd) %>% 
    mutate(iter = 1:1000)
  fxd
})
names(sim2_scenario_res) <- str_sub(sim2_scenarios, -46, -5)
sim2_scenario_res <- bind_rows(sim2_scenario_res, .id = "scenario")



# Find the iteration when got the mean value (or closest to it)
mean_scenario <- scenario_res %>% 
  group_by(scenario) %>% 
  mutate(min_diff = abs(mean - mean(mean))) %>% 
  arrange(scenario, min_diff) %>% 
  slice(1) %>% 
  ungroup()

# Take the mean value and qintiles
scenario_res_q <- scenario_res %>% 
  select(scenario, mean, `0.025quant`, `0.975quant`) %>% 
  group_by(scenario) %>% 
  summarise_all(.funs = list(est = "mean",
                          lci =function(x) quantile(x, 0.025),
                          uci =function(x) quantile(x, 0.975))) %>% 
  ungroup()

scenario_res_q <- scenario_res_q %>% 
  gather(key = "scenario_new", value = "value", -scenario) %>% 
  separate(scenario_new, into = c("stat", "smry"), sep = "_")

## Separate names
scenario_names_lng <- ScenarioNames(scenario_res_q$scenario)
scenario_res_q <- bind_cols(scenario_names_lng, scenario_res_q)

## Spreadstatistic to wide
scenario_res_q <- scenario_res_q %>% 
  spread(key = smry, value = value)

scenario_res2 <- scenario_res_q %>% 
  mutate(trial_drug = trial + atc5,
         result = factor(stat, levels = c("0.025quant", "mean", "0.975quant")),
         my_alpha = if_else(result == "mean", 1, 0),
         trial = paste0("Trial ", trial),
         drug = paste0("Drug ", drug),
         atc5 = paste0("Class ", atc5)) %>% 
  as_tibble()

pd <- position_dodge(width = 0.4)

emphasise_class <- ggplot(scenario_res2,
                aes(x = atc5, y = est, ymin = lci, ymax = uci, colour = result,
                    alpha = my_alpha)) +
  geom_errorbar(position = pd) +
  geom_point(position = pd) +
  facet_grid(trial ~ drug) +
  scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
  scale_y_continuous("Effect estimate") +
  scale_alpha(range = c(0.4, 1), guide = FALSE) +
  scale_colour_discrete("") 

emphasise_trial <- emphasise_class %+% scenario_res2 +
  aes(x = trial) +
  facet_grid(atc5 ~ drug)

emphasise_drug <- emphasise_class %+% scenario_res2 +
  aes(x = drug) +
  facet_grid(atc5 ~ trial)

pdf("figures/unix_sim1.pdf")
emphasise_class + ggtitle("Drug class variation on x-axis")
emphasise_drug + ggtitle("Drug variation on x-axis")
emphasise_trial + ggtitle("Trial variation on x-axis")
dev.off()

# all_single_plot <- ggplot(scenario_res2,
#                 aes(x = atc5, y = q50, ymin = q2.5, ymax = q97.5, colour = trial_drug,
#                     alpha = my_alpha, group = result)) +
#   geom_errorbar(position = pd) +
#   geom_point(position = pd) +
#   scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
#   scale_y_continuous("Effect estimate") +
#   scale_alpha(range = c(0.4, 1), guide = FALSE)

# Legacy code?

tiff("figures/variation_estimates.tiff", res = 600, compression = "lzw",
     unit = "in", height = 7, width = 7)
emphasise_class + scale_x_discrete("Between class variation", labels = c(0.05, 0.15, 0.25))
dev.off()


## Examine predictive distributions
mean_scenario <- bind_cols(scenarios_names, mean_scenario)
saveRDS(mean_scenario, "scratch_data/mean_scenario")
## RE-run models for 05_05_05 and 
