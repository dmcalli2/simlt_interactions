# Scratch exploration of plotting scenario results

library(tidyverse)
library(stringr)
library(ggplot2)


scenario_res <- readRDS("scratch_data/scenario_res.Rds")

# Ascertain total variation for each scenario and sort by this

scen_res_sum <- scenario_res %>%
  group_by(scenario) %>%
  summarise(total_vrn = mean(sd)) %>%
  inner_join(scenario_res %>% 
               select(scenario,sim,como_prev)%>% 
               distinct()) %>%
  group_by(sim,como_prev) %>%
  arrange(total_vrn ,.by_group=TRUE)%>%
  ungroup()

# Ascertain minimum, maximum, and median variation scenario for each sim

key_scens_med <- scen_res_sum %>%
  group_by(sim, como_prev)%>%
  summarise(total_vrn = median(total_vrn)) %>%
  ungroup() %>%
  inner_join(scen_res_sum %>% select(scenario,total_vrn, como_prev)) %>%
  mutate(key_scen_type = "med_total_vrn")
key_scens_min <- scen_res_sum %>%
  group_by(sim, como_prev)%>%
  slice(which.min(total_vrn))%>%
  mutate(key_scen_type = "min_total_vrn")%>%
  ungroup() 
key_scens_max <- scen_res_sum %>%
  group_by(sim, como_prev)%>%
  slice(which.max(total_vrn))%>%
  mutate(key_scen_type = "max_total_vrn")%>%
  ungroup() 

key_scens <- key_scens_max %>%
  bind_rows(key_scens_med) %>%
  bind_rows(key_scens_min)

saveRDS(key_scens, "scratch_data/key_scens")

scen_res_sum <- scen_res_sum %>%
  left_join(key_scens) %>%
  filter(como_prev=="std")

#Visualise

ggplot(scen_res_sum,
       aes(x = reorder(scenario,total_vrn), y = total_vrn, fill = key_scen_type)) +
  geom_bar(stat="identity") +
  facet_wrap(~sim, scales="free") +
  coord_flip()

# Explore pirate plots to show variation across iterations within (key) scenarios for std prevalence

key_scens_std <- key_scens %>%
  filter(como_prev == "std")

key_scens_list <- as.vector(key_scens_std$scenario)

key_scens_all_its <- scenario_res %>%
  filter(como_prev == "std", scenario %in% key_scens_list ) %>%
  left_join(key_scens_std)

# Take the mean value and qintiles
scenario_res_q <- key_scens_all_its %>% 
  select(scenario, mean, `0.025quant`, `0.975quant`) %>% 
  group_by(scenario) %>% 
  summarise_all(.funs = list(est = "mean",
                             lci =function(x) quantile(x, 0.025),
                             uci =function(x) quantile(x, 0.975))) %>% 
  ungroup() %>%
  inner_join(key_scens_std %>% select(scenario,sim)) %>%
  distinct()


ggplot()+
  geom_jitter(data=key_scens_all_its, aes(x = scenario, y = mean, color = key_scen_type)) +
  geom_violin(data=key_scens_all_its, aes(x = scenario, y = mean ),fill=NA, size=1 ) +
  geom_point(data=scenario_res_q, aes(x = scenario, y = mean_est, size = 1)) +
  geom_errorbar(data=scenario_res_q,  aes(x = scenario, ymin = mean_lci, ymax = mean_uci, width=0.3 )) +
  facet_wrap(~sim, scales="free_x") 

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

# For same select scenarios, show effect of como_prev

key_scens_all_como_all_its <- key_scens %>%
  left_join(scenario_res)
# Take the mean value and qintiles
scenario_res_q_como <- key_scens_all_como_all_its %>% 
  select(scenario,como_prev, mean, `0.025quant`, `0.975quant`) %>% 
  group_by(scenario, como_prev) %>% 
  summarise_all(.funs = list(est = "mean",
                             lci =function(x) quantile(x, 0.025),
                             uci =function(x) quantile(x, 0.975))) %>% 
  ungroup() %>%
  inner_join(key_scens %>% select(scenario,sim,como_prev,key_scen_type)) %>%
  distinct()

ggplot()+
  geom_jitter(data=key_scens_all_como_all_its, aes(x = key_scen_type, y = mean, color = key_scen_type)) +
  geom_violin(data=key_scens_all_como_all_its, aes(x = key_scen_type, y = mean ),fill=NA, size=1 ) +
  geom_point(data=scenario_res_q_como, aes(x = key_scen_type, y = mean_est, size = 1)) +
  geom_errorbar(data=scenario_res_q_como,  aes(x = key_scen_type, ymin = mean_lci, ymax = mean_uci, width=0.3 )) +
  facet_grid(como_prev~sim, scales="free") 

# Same idea, but also plotting estimates of quantiles as well

key_scens_all_como_all_its_lng <- key_scens_all_como_all_its %>%
  rename(loq = `0.025quant`, upq = `0.975quant`) %>%
  select(scenario, mean, loq, upq) %>%
  gather(key = "val_type", value = "value", -scenario)

key_scens_all_como_all_its_lng <- key_scens_all_como_all_its_lng%>%
  inner_join(key_scens %>% select(scenario, sim,como_prev, key_scen_type) %>% distinct())%>%
  filter(como_prev %in% c("hi","lo"))


scenario_res_q_como2 <-  scenario_res_q_como %>%
  gather(key = "scenario_new", value = "value", -scenario,-como_prev,-sim,-key_scen_type) %>% 
  separate(scenario_new, into = c("stat", "smry"), sep = "_") %>% 
  spread(key = smry, value = value) %>%
  select(-key_scen_type) %>%
  inner_join(key_scens %>% select(scenario, sim,como_prev, key_scen_type) %>% distinct()) %>%
  filter(como_prev %in% c("hi","lo"))


key_scens_all_como_all_its_lng$key_scen_type <- as.factor(key_scens_all_como_all_its_lng$key_scen_type )
key_scens_all_como_all_its_lng$key_scen_type <- factor(key_scens_all_como_all_its_lng$key_scen_type,levels(key_scens_all_its_lng$key_scen_type)[c(3,2,1)])                  
scenario_res_q_como2$key_scen_type <- as.factor(scenario_res_q_como2$key_scen_type )
scenario_res_q_como2$key_scen_type <- factor(scenario_res_q_como2$key_scen_type,levels(scenario_res_q2$key_scen_type)[c(3,2,1)])                  

key_scens_all_como_all_its <- key_scens_all_como_all_its %>%
  filter(como_prev %in% c("hi","lo"))

scen_res_sum$sim <- as.factor(scen_res_sum$sim)
levels(scen_res_sum$sim)[levels(scen_res_sum$sim)==1] <- "Diabetes"
levels(scen_res_sum$sim)[levels(scen_res_sum$sim)==2]   <- "Inflammatory conditions"
# Show effect (or lack of effect) of comorbidity prevalence
ggplot()+
  geom_jitter(data=key_scens_all_como_all_its_lng, aes(x = key_scen_type, y = value, color = val_type), width = 0.4 )+
  geom_errorbar(data=scenario_res_q_como2,  aes(x = key_scen_type, ymin = lci, ymax = uci), width=0.1, size=1, color="grey20", position=position_dodge(0.001) ) +
  geom_violin(data=key_scens_all_como_all_its, aes(x = key_scen_type, y = mean),alpha=0.3, size=1.2, colour = "springgreen4" )+
  geom_point(data=scenario_res_q_como2, aes(x = key_scen_type, y = est), size=2, shape=2, fill="grey20") +
  facet_wrap(~sim+como_prev, scales="free_y") + 
  theme(axis.text.x = element_text(angle = 0), 
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 20, r = 0, b = 20, l = 0)),
        text =element_text(size = 14.5),
        panel.background = element_rect(fill = "white", colour = "grey80"))
