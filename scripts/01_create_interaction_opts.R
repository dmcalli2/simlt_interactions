# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next Calculates overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)
library(tidyverse)
## Read in diabetes trials
load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")

# Remove single alpha-glucosidase inhibitor trial
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")


## Make same answer each set of classes (not each class)
set.seed(1234)

### Each simulation scenario, overall effects
main_scen <- expand.grid(
  cept = 0,
  como = c(-0.2,-0.1, 0, 0.1, 0.2),
  allc = c(-0.2,-0.1, 0, 0.1, 0.2),
  actn = c(-0.2,-0.1, 0, 0.1, 0.2)
)
### Each simulation scenario, variation around the overall effect
varn_scen <- list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25))

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all drugs, classes and trial

# Simulate variation around wider group (A10B) effect at trial, drug and class level
SampleVarn <- function(item, smpl_lvls, n_iter = 1000){
  my_mtrx <- matrix(nrow = length(item) * n_iter,
                    ncol = length(smpl_lvls))
  
  for(j in seq_along(smpl_lvls)){
    my_mtrx[,j] <- rnorm(length(item) * n_iter, 0, smpl_lvls[j])
  }
  rownames(my_mtrx) <- rep(item, n_iter)
  colnames(my_mtrx) <- smpl_lvls
  my_mtrx
}

trial <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(diabetes_final$nct_id, x)
                }
)

drug <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$drug), x)
                }
)

atc5 <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(diabetes_final$atc_5), x)
                }
)

## Combine variation around A10B effects for intercept, comorbidity, allocation and itneraction
## into single estimate from trial, drug and class (A10BA, A10BB etc)
atc5 <- atc5$actn
atc5_names <- rownames(atc5)
colnames(atc5) <- paste("atc5", colnames(atc5), sep = "_")
atc5 <- as_tibble(atc5) %>%
mutate(atc_5 = atc5_names)

drug <- drug$actn
drug_names <- rownames(drug)
colnames(drug) <- paste("drug", colnames(drug), sep = "_")
drug <- as_tibble(drug) %>%
mutate(drug = drug_names)

trial <- trial$actn
trial_names <- rownames(trial)
colnames(trial) <- paste("trial", colnames(trial), sep = "_")
trial <- as_tibble(trial) %>%
mutate(nct_id = trial_names)


trial <- trial %>% 
  arrange(nct_id)

drug <- drug %>% 
  inner_join(diabetes_final %>%  select(nct_id, drug)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

atc5 <- atc5 %>% 
  inner_join(diabetes_final %>%  select(nct_id, atc_5)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

## Combine all interaction effects
actn <- bind_cols(atc5, trial, drug) %>% 
  arrange(atc_5, drug, nct_id) %>% 
  mutate(iteration = rep(1:1000, 161)) %>% 
    select(atc_5, drug, nct_id, iteration, everything())

## Examine interaction effects at drug-class level
actn %>%
  select(iteration, atc_5, starts_with("atc5")) %>% 
  distinct() %>%
  # Examine SD across classes for each iteration
  group_by(iteration) %>% 
  summarise_at(vars(starts_with("atc5")), sd) %>%
  # Summarise results across all iterations
  summarise_at(vars(starts_with("atc5")), mean)

saveRDS(actn, file = "scratch_data/interactn_opts.Rds")
save(diabetes_final, file = "scratch_Data/interaction_opts_ordering.Rds")
