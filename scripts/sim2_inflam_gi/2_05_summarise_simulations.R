# 03_summarise_simulations
library(tidyverse)
library(stringr)
library(ggplot2)

# Identify saved results
scenarios <- list.files("unix_results/sim2/sd1/", patt = "scen")

ScenarioNames <- function (scenarios) {
  scenarios_names <- str_split_fixed( str_sub(scenarios, 1, -5), "_", n = 10)
  scenarios_names <- apply(scenarios_names[, c(4, 6, 8, 10)], 2, as.double)
  scenarios_names <- as.data.frame(scenarios_names)
  names(scenarios_names) <- c("path", "moa", "trial", "drug")
  scenarios_names
}

scenarios_names <- ScenarioNames(scenarios)

# read and convert to data frame for each scenario
scenario_res <- lapply(scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("unix_results/sim2/sd1/", each_scenario))})

scenario_res <- map(scenario_res, function (scen){
  fxd <- map(scen, ~ .x$fixed) 
  fxd <- do.call(rbind, fxd)
  fxd <- as.tibble(fxd) %>% 
  mutate(iter = 1:1000)
  fxd
})
names(scenario_res) <- scenarios
scenario_res <- bind_rows(scenario_res, .id = "scenario")

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
  mutate(trial_drug = trial + drug,
         result = factor(stat, levels = c("0.025quant", "mean", "0.975quant")),
         my_alpha = if_else(result == "mean", 1, 0),
         trial = paste0("Trial ", trial),
         drug = paste0("Drug ", drug),
         moa = paste0("MoA ", moa),
         path = paste0("Broad_pth ", path)) %>% 
  as_tibble()

pd <- position_dodge(width = 0.8)

  emphasise_class <- ggplot(scenario_res2,
                aes(x = moa, y = est, ymin = lci, ymax = uci, colour = result, shape = path, 
                    alpha = my_alpha)) +
  geom_errorbar(position = pd) +
  geom_point(position = pd) +
  facet_grid(trial ~ drug) +
  scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
  scale_y_continuous("Effect estimate") +
  scale_alpha(range = c(0.4, 1), guide = FALSE) +
  theme_light() + 
  scale_colour_discrete("") 
  emphasise_class
  
emphasise_trial <- emphasise_class %+% scenario_res2 +
  aes(x = trial) +
  facet_grid(moa ~ drug)

emphasise_drug <- emphasise_class %+% scenario_res2 +
  aes(x = drug) +
  facet_grid(moa ~ trial)

pdf("figures/unix_sim2_sd1.pdf")
emphasise_class + ggtitle("MoA variation on x-axis")
emphasise_drug + ggtitle("Drug variation on x-axis")
emphasise_trial + ggtitle("Trial variation on x-axis")
dev.off()

# all_single_plot <- ggplot(scenario_res2,
#                 aes(x = moa, y = q50, ymin = q2.5, ymax = q97.5, colour = trial_drug,
#                     alpha = my_alpha, group = result)) +
#   geom_errorbar(position = pd) +
#   geom_point(position = pd) +
#   scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
#   scale_y_continuous("Effect estimate") +
#   scale_alpha(range = c(0.4, 1), guide = FALSE)

## Legacy code?

tiff("figures/variation_estimates_sim2_sd1.tiff", res = 600, compression = "lzw",
     unit = "in", height = 7, width = 7)
emphasise_class + scale_x_discrete("Between class variation", labels = c(0.05, 0.15, 0.25))
dev.off()


## Examine predictive distributions
mean_scenario <- bind_cols(scenarios_names, mean_scenario)
saveRDS(mean_scenario, "scratch_data/mean_scenario_sim2_sd1")
## RE-run models for 05_05_05 and 
