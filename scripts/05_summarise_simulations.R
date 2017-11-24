# 03_summarise_simulations
library(tidyverse)
library(stringr)
library(ggplot2)

# scenarios <- c("scenario_atc5_0.05_trial_0.05_drug_0.05.rds",
#                "scenario_atc5_0.1_trial_0.1_drug_0.1.rds",
#                "scenario_atc5_0.1_trial_0.1_drug_0.25.rds",
#                "scenario_atc5_0.1_trial_0.25_drug_0.1.rds", 
#                "scenario_atc5_0.25_trial_0.1_drug_0.1.rds",
#                "scenario_atc5_0.25_trial_0.25_drug_0.25.rds")


scen1 <- readRDS("unix_results/scenario_atc5_0.05_trial_0.05_drug_0.05.rds")
scen2 <- readRDS("unix_results/scenario_atc5_0.25_trial_0.25_drug_0.25.rds")

scenarios <- list.files("unix_results/", patt = "scen")

scenarios_names <- str_split_fixed( str_sub(scenarios, 1, -5), "_", n = 7)
scenarios_names <- apply(scenarios_names[, c(3, 5, 7)], 2, as.double)
scenarios_names <- as.data.frame(scenarios_names)
names(scenarios_names) <- c("atc5", "trial", "drug")

scenario_res <- lapply(scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("unix_results/", each_scenario))
  fxd <- lapply(each_scenario, function(x) x$fixed)
  fxd <- do.call(rbind, fxd)
  fxd <- as.data.frame(fxd)
  sapply(fxd, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
})
names(scenario_res) <- scenarios

scenario_res <- lapply(seq_along(scenarios), function(i) {
  scenario_res <- as.data.frame(scenario_res[[i]])
  scenario_res <- cbind(scenario_res, scenarios_names[i,])
  scenario_res$qs <- paste0("q", c(2.5, 50, 97.5))
  scenario_res
})

scenario_res <- bind_rows(scenario_res)


scenario_res2 <- scenario_res %>% 
  gather(key = "var", value = "value", -atc5, -drug, -trial, -qs) %>% 
  spread(key = qs, value = value) %>% 
  filter(var %in% c("0.025quant", "mean", "0.975quant")) %>% 
  mutate(result = factor(var, c("0.025quant", "mean", "0.975quant"),
                         c("lci", "mean", "uci")),
         my_alpha = if_else(result == "mean", 1, 0),
         trial = paste0("Trial ", trial),
         drug = paste0("Drug ", drug),
         atc5 = paste0("Class ", atc5))

pd <- position_dodge(width = 0.4)

emphasise_class <- ggplot(scenario_res2,
                aes(x = atc5, y = q50, ymin = q2.5, ymax = q97.5, colour = result,
                    alpha = my_alpha)) +
  geom_errorbar(position = pd) +
  geom_point(position = pd) +
  facet_grid(trial ~ drug) +
  scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
  scale_y_continuous("Effect estimate") +
  scale_alpha(range = c(0.4, 1), guide = FALSE)

emphasise_trial <- emphasise_class %+% scenario_res2 +
  aes(x = trial) +
  facet_grid(atc5 ~ drug)

emphasise_drug <- emphasise_class %+% scenario_res2 +
  aes(x = drug) +
  facet_grid(atc5 ~ trial)

pdf("unix_sim1.pdf")
emphasise_class + ggtitle("Drug class variation on x-axis")
emphasise_drug + ggtitle("Drug variation on x-axis")
emphasise_trial + ggtitle("Trial variation on x-axis")
dev.off()
