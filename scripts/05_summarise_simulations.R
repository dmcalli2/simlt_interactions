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
  names(scenarios_names) <- c("como_prev","atc5_moa", "trial", "drug")
  }
  else if(sim == 2) {scenarios_names <-  str_split_fixed( str_sub(scenarios, -46, -5), "_", n = 9)
  scenarios_names <- cbind(scenarios_names[,1],apply(scenarios_names[, c(3, 5, 7, 9)], 2, as.double))
  scenarios_names <- as.data.frame(scenarios_names)
  names(scenarios_names) <- c("como_prev","path","atc5_moa", "trial", "drug")
  }
  scenarios_names
}

sim1_scenarios_names <- ScenarioNames(sim1_scenarios,sim=1)
sim2_scenarios_names <- ScenarioNames(sim2_scenarios,sim=2)

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

#Tidy up and combine

sim1_scenario_res <- sim1_scenario_res %>%
  mutate(como_prev = str_sub(scenario, 1, 2),
         atc_moa = str_sub(scenario, 9, 12),
         trial = str_sub(scenario, 20,23),
         drug = str_sub(scenario, 30, 33),
         path = "0",
         sim = 1)

sim2_scenario_res <- sim2_scenario_res %>%
  mutate(como_prev = str_sub(scenario, 1, 2),
         path = str_sub(scenario, 9, 12),
         atc_moa = str_sub(scenario, 18,21),
         trial = str_sub(scenario, 29, 32),
         drug = str_sub(scenario, 39, 42),
         sim = 2)

scenario_res <- sim1_scenario_res %>%
  bind_rows(sim2_scenario_res) %>%
  mutate(como_prev = ifelse(como_prev %in% "td", "std", como_prev)) %>%
  mutate(scenario = paste0("cp",scenario)  )

saveRDS(scenario_res, "scratch_data/scenario_res.Rds")

##### Review continous como_prev results

como_prev_scen_res <- scenario_res %>%
  filter(!como_prev %in% c("hi","lo","std")) %>%
  mutate(como_prev = ifelse(!como_prev %in% c(".2",".1"), paste0(".",como_prev), como_prev ))
como_prev_scen_res$como_prev <- as.numeric(como_prev_scen_res$como_prev)

# Find the iteration when got the mean value (or closest to it)
mean_scenario <- como_prev_scen_res %>% 
  group_by(scenario) %>% 
  mutate(min_diff = abs(mean - mean(mean))) %>% 
  arrange(scenario, min_diff) %>% 
  slice(1) %>% 
  ungroup()

# Take the mean value and qintiles
scenario_res_q <- como_prev_scen_res %>% 
  select(scenario, mean, `0.025quant`, `0.975quant`) %>% 
  group_by(scenario) %>% 
  summarise_all(.funs = list(est = "mean",
                             lci =function(x) quantile(x, 0.025),
                             uci =function(x) quantile(x, 0.975))) %>% 
  ungroup()

scenario_res_q <- scenario_res_q %>% 
  gather(key = "scenario_new", value = "value", -scenario) %>% 
  separate(scenario_new, into = c("stat", "smry"), sep = "_")%>%
  left_join(como_prev_scen_res %>%
              distinct(sim, como_prev, path, atc_moa, drug, trial, .keep_all=TRUE) %>%
              select(scenario,sim, como_prev, path, atc_moa, drug, trial))

## Spreadstatistic to wide
scenario_res_q <- scenario_res_q %>% 
  spread(key = smry, value = value)

scenario_res2 <- scenario_res_q %>% 
  mutate(result = factor(stat, levels = c("0.025quant", "mean", "0.975quant")),
         my_alpha = if_else(result == "mean", 1, 0),
         Trial = paste0("Trial ", trial),
         Drug = paste0("Drug ", drug),
         Atc_moa = paste0("Class ", atc_moa),
         path = as.numeric(path),
         atc_moa = as.numeric(atc_moa),
         drug = as.numeric(drug),
         trial = as.numeric(trial),
         scenario = str_sub(scenario, 6,45)) %>%    
        as_tibble()                     

pd <- position_dodge(width = 0)

scenario_res2$path <- as.character(scenario_res2$path)

labels <- c(atc5_0.05_trial_0.05_drug_0.05 = "Diabetes: \nLow variation", atc5_0.25_trial_0.25_drug_0.25= "Diabetes: \nHigh variation", 
            path_0.05_moa_0.05_trial_0.05_drug_0.05= "Rheumatoid: \nLow variation", path_0.25_moa_0.25_trial_0.25_drug_0.25 = "Rheumatoid: \nHigh variation")

emphasise_class <- ggplot(scenario_res2,
                          aes(x = como_prev, y = est, ymin = lci, ymax = uci, colour = result,
                              alpha = my_alpha)) +
  geom_errorbar(position = pd) +
  geom_point(position = pd) +
  facet_grid(.~  scenario, labeller=labeller(scenario = labels) ) +
  scale_alpha(range = c(0.4, 1), guide = FALSE) + 
  scale_colour_discrete("Value type \n(average across all \niterations of a scenario, \ndrawn from posterior \ndistribution of \neffect estimate):",
                        breaks=c("0.025quant","mean","0.975quant"),
                        labels=c("2.5th percentile","Mean","97.5th percentile") )+
  geom_hline( aes(yintercept=-0.1, linetype = "Simulated effect \nat WDG level"), color = "black")+
  theme(axis.text.x = element_text(angle=0), 
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        strip.text = element_text(size=11),
        legend.title = element_text(size =11)) +
  scale_y_continuous("Interaction effect estimate",breaks=seq(-1.5,1.5,0.5)) +
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c( "black")))) + 
  scale_x_continuous("Simulated comorbidity prevalence",breaks=seq(0,0.2,0.05), labels = c('','0.05','','0.15','')) 

emphasise_class
tiff("figures/effect_como_prev.tiff", res = 600, compression = "lzw", unit = "in",
     height = 5, width = 8)
emphasise_class
dev.off()

#### Review std como_prev results

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
scenario_names_lng <- ScenarioNames(scenario_res_q$scenario, sim=1)
scenario_res_q <- scenario_res_q %>%
  left_join(scenario_res %>%
              distinct(sim, como_prev, path, atc_moa, drug, trial, .keep_all=TRUE) %>%
              select(scenario,sim, como_prev, path, atc_moa, drug, trial))

scenario_res_q$como_prev <- as.factor(scenario_res_q$como_prev)


## Spreadstatistic to wide
scenario_res_q <- scenario_res_q %>% 
  spread(key = smry, value = value)

scenario_res2 <- scenario_res_q %>% 
  mutate(result = factor(stat, levels = c("0.025quant", "mean", "0.975quant")),
         my_alpha = if_else(result == "mean", 1, 0),
         Trial = paste0("Trial ", trial),
         Drug = paste0("Drug ", drug),
         Atc_moa = paste0("Class ", atc_moa),
         path = as.numeric(path),
         atc_moa = as.numeric(atc_moa),
         drug = as.numeric(drug),
         trial = as.numeric(trial)) %>% 
  filter(como_prev %in% c("std")) %>%    ##### Comorbidity prevalence appears to make little difference
  as_tibble()                     #####  so drop hi/lo here and use continuous como_prev to investigate separately

pd <- position_dodge(width = 1)

scenario_res2$path <- as.character(scenario_res2$path)

emphasise_class <- ggplot(scenario_res2,
                          aes(x = Atc_moa, y = est, ymin = lci, ymax = uci, colour = result,
                              alpha = my_alpha, shape = como_prev)) +
  geom_errorbar(position = pd) +
  geom_point(position = pd) +
  facet_grid(sim+ path ~  Trial + Drug  ) +
  scale_x_discrete("", labels = c(0.05, 0.15, 0.25)) +
  scale_y_continuous("Effect estimate") +
  scale_alpha(range = c(0.4, 1), guide = FALSE) + 
  scale_colour_discrete("") 

emphasise_class


# Ascertain total variation for each scenario and sort by this

scen_res_sum <- scenario_res %>%
  group_by(scenario) %>%
  summarise(total_vrn = mean(sd)) %>%
  inner_join(scenario_res_q %>% 
               select(scenario,sim,como_prev,est,stat,lci,uci)) %>%
  filter(como_prev=="std")

scen_res_sum$sim <- as.factor(scen_res_sum$sim)
levels(scen_res_sum$sim)[levels(scen_res_sum$sim)==1] <- "Diabetes trial-set"
levels(scen_res_sum$sim)[levels(scen_res_sum$sim)==2]   <- "Rheumatoid trial-set"
  
show_all <- ggplot() +
  geom_point(data=scen_res_sum, aes(x=reorder(scenario,total_vrn),y=est, colour=stat))+
  geom_errorbar(data=scen_res_sum, aes(x= reorder(scenario,total_vrn), ymax=uci,ymin=lci, colour= stat)) +
  facet_wrap(~sim, scales= "free_x")  +
  scale_colour_discrete("Value type \n(average across all \niterations of a scenario, \ndrawn from posterior \ndistribution of \neffect estimate):",
                        breaks=c("0.025quant","mean","0.975quant"),
                        labels=c("2.5th percentile","Mean","97.5th percentile") )+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"),
        panel.grid.major = element_line(colour = "grey80"),
        strip.text = element_text(size=11),
        legend.title = element_text(size =11)) +
  scale_y_continuous(breaks=seq(-1.6,1.6,0.2)) +
  scale_x_discrete(breaks = NULL)+
  xlab("Simulated scenarios, in order of \nincreasing total variation in the hierarchy")+
  ylab("Interaction effect estimate")+ 
  geom_hline(data=scen_res_sum, aes(yintercept=-0.1, linetype = "Simulated effect \nat WDG level"), color = "black")+
  scale_linetype_manual(name = "", values = c(2), 
                        guide = guide_legend(override.aes = list(color = c( "black")))) 

show_all
tiff("figures/recovered_effects_main.tiff", res = 600, compression = "lzw", unit = "in",
     height = 8, width = 8)
show_all
dev.off()

## Tablulate precision of effect est recovery for key scens 

scenarios <- c('cptd_atc5_0.05_trial_0.05_drug_0.05','cptd_atc5_0.15_trial_0.15_drug_0.15','cptd_atc5_0.25_trial_0.25_drug_0.25',
               'cptd_path_0.05_moa_0.05_trial_0.05_drug_0.05','cptd_path_0.15_moa_0.15_trial_0.15_drug_0.15','cptd_path_0.25_moa_0.25_trial_0.25_drug_0.25')

scen_res_sum_key <- scen_res_sum %>%
  filter(scenario %in% scenarios) %>%
  select(scenario,sim,est,stat) %>%
  spread(key = stat, value = est) 

write.csv(scen_res_sum_key, file = "./tables/key_scens_sum.csv")


## Show comparison with model excluding DC levels

sim_scenarios <- list.files("unix_results/noclasscomp/", recursive=TRUE, pattern = "")

# read and convert to data frame for each scenario
nodc_scenario_res <- lapply(sim_scenarios, function(each_scenario){
  each_scenario <- readRDS(paste0("unix_results/noclasscomp/", each_scenario))})

nodc_scenario_res <- map(nodc_scenario_res, function (scen){
  fxd <- map(scen, ~ .x$fixed) 
  fxd <- do.call(rbind, fxd)
  fxd <- as.tibble(fxd) %>% 
    mutate(iter = 1:1000)
  fxd
})
names(nodc_scenario_res) <- str_sub(sim_scenarios, 11, -5)
nodc_scenario_res <- bind_rows(nodc_scenario_res, .id = "scenario") 
nodc_scenario_res <- nodc_scenario_res %>%
  mutate(model = "NoDC")
scens <-  unique(nodc_scenario_res$scenario)

scenario_res <- readRDS("scratch_data/scenario_res") 

scen_res_both <- nodc_scenario_res %>%
  bind_rows(scenario_res %>% 
              filter(scenario %in% scens) %>%
              mutate(model = "Full") %>%
              select(scenario, mean,sd,`0.025quant`,`0.5quant`,`0.975quant`,mode,kld,iter,model )) %>%
  mutate(sim = ifelse(scenario %in% scens[seq(1,(length(scens)/2),1)], 1,2))


ggplot()+
  geom_jitter(data=scen_res_both, aes(x = scenario, y = mean, color = model)) +
  geom_violin(data=scen_res_both, aes(x = scenario, y = mean ),fill=NA, size=1 ) +
  facet_grid(model~sim, scales="free_x") 

# Same idea, but also plotting estimates of quantiles as well

key_scens_all_its_lng <- key_scens_all_its %>%
  rename(loq = `0.025quant`, upq = `0.975quant`) %>%
  select(scenario, mean, loq, upq) %>%
  gather(key = "val_type", value = "value", -scenario)

key_scens_all_its_lng <- key_scens_all_its_lng%>%
  inner_join(key_scens_std %>% select(scenario, sim, key_scen_type) %>% distinct())

scenario_res_q2 <-  scenario_res_q %>%
  gather(key = "scenario_new", value = "value", -scenario) %>% 
  separate(scenario_new, into = c("stat", "smry"), sep = "_") %>% 
  spread(key = smry, value = value) %>%
  filter(stat != "sim") %>%
  inner_join(key_scens_std %>% select(scenario, sim, key_scen_type )) %>%
  distinct()

key_scens_all_its_lng$key_scen_type <- as.factor(key_scens_all_its_lng$key_scen_type )
key_scens_all_its_lng$key_scen_type <- factor(key_scens_all_its_lng$key_scen_type,levels(key_scens_all_its_lng$key_scen_type)[c(3,2,1)])                  
scenario_res_q2$key_scen_type <- as.factor(scenario_res_q2$key_scen_type )
scenario_res_q2$key_scen_type <- factor(scenario_res_q2$key_scen_type,levels(scenario_res_q2$key_scen_type)[c(3,2,1)])                  
key_scens_all_its_lng$sim <- as.factor(key_scens_all_its_lng$sim)
levels(key_scens_all_its_lng$sim)[levels(key_scens_all_its_lng$sim)==1] <- "Diabetes"
levels(key_scens_all_its_lng$sim)[levels(key_scens_all_its_lng$sim)==2]   <- "Inflammatory conditions"
scenario_res_q2$sim <- as.factor(scenario_res_q2$sim)
levels(scenario_res_q2$sim)[levels(scenario_res_q2$sim)==1] <- "Diabetes"
levels(scenario_res_q2$sim)[levels(scenario_res_q2$sim)==2]   <- "Inflammatory conditions"
key_scens_all_its$sim <- as.factor(key_scens_all_its$sim)
levels(key_scens_all_its$sim)[levels(key_scens_all_its$sim)==1] <- "Diabetes"
levels(key_scens_all_its$sim)[levels(key_scens_all_its$sim)==2]   <- "Inflammatory conditions"

ggplot()+
  geom_jitter(data=key_scens_all_its_lng, aes(x = key_scen_type, y = value, color = val_type), width = 0.4 )+
  geom_errorbar(data=scenario_res_q2,  aes(x = key_scen_type, ymin = lci, ymax = uci), width=0.1, size=1, color="grey20", position=position_dodge(0.001) ) +
  geom_violin(data=key_scens_all_its, aes(x = key_scen_type, y = mean),alpha=0.3, size=1.2, colour = "springgreen4" )+
  geom_point(data=scenario_res_q2, aes(x = key_scen_type, y = est), size=2, shape=2, fill="grey20") +
  facet_wrap(~sim) +
  theme( axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
         axis.title.y = element_text(margin = margin(t = 0, r = 30, b = , l = 0)),
         text =element_text(size = 14.5),
         panel.background = element_rect(fill = "white", colour = "grey80")) +
  xlab("Selected scenarios, all iterations")+
  ylab("Interaction effect estimate from model")+ 
  scale_x_discrete(labels= c("Minimum variation\n scenario","Median variation\n scenario","Maximum variation\n scenario" ))+
  geom_hline(yintercept = -0.1, linetype = "dashed", color = "black")



names(nodc_scenario_res)## Examine predictive distributions

saveRDS(mean_scenario, "scratch_data/mean_scenario")

########## modified for both sims to here

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
