# 11_simulate_data_base_r
# Initially generates random samples from normal distribution. THese will be the same for each scenario
# Next Calculates overall group effect and variance for each component that will make-up the final scenarios
# eg (intercept, comorbidity effect, treatment allocation and interaction)
library(tidyverse)
## Read in inflam/GI trials
load("data/rheum_final.Rdata")

table(rheum_final$conditions, rheum_final$atc_5)

## Make same answer each set of classes (not each class)
set.seed(1234)


### Each simulation scenario, variation around the overall effect
varn_scen <- list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25))

### Determine number of samples needed for each scenario, do random sampling once for all analyses
### for all trials, drugs, MoAs, and broader pathways  ### To decide- what to do about indication?

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
                 SampleVarn(unique(rheum_final$nct_id), x)
               }
)

drug <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(rheum_final$drug), x)
                }
)

moa <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(rheum_final$moa), x)
                }
)

path <- lapply(list(cept = c(0.25, 0.5),
                     como = c(0.25, 0.5),
                     allc = c(0.05, 0.10, 0.15, 0.20, 0.25),
                     actn = c(0.05, 0.10, 0.15, 0.20, 0.25)),
                function(x) {
                  SampleVarn(unique(rheum_final$brd_drug_pth), x)
                }
)
## Combine variation around overall effects for intercept, comorbidity, allocation and interaction 
## into single estimate from trial, drug and class (A10BA, A10BB etc)
path <- path$actn
path_names <- rownames(path)
colnames(path) <- paste("path", colnames(path), sep = "_")
path <- as_tibble(path) %>%
mutate(brd_drug_pth = path_names)

moa <- moa$actn
moa_names <- rownames(moa)
colnames(moa) <- paste("moa", colnames(moa), sep = "_")
moa <- as_tibble(moa) %>%
  mutate(moa = moa_names)

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
  inner_join(rheum_final %>%  select(nct_id, drug, conditions)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

moa <- moa %>% 
  inner_join(rheum_final %>%  select(nct_id, moa)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

path <- path %>% 
  inner_join(rheum_final %>%  select(nct_id, brd_drug_pth)) %>% 
  arrange(nct_id) %>% 
  select(-nct_id)

## Combine all interaction effects
actn <- bind_cols(path, moa, drug, trial) %>% 
  arrange(brd_drug_pth, moa, drug,  nct_id) %>% 
  mutate(iteration = rep(1:1000, length(rheum_final$nct_id))) %>% 
    select(brd_drug_pth, moa, drug,indication = conditions, iteration, everything())


## Examine interaction effects at drug-class (MoA) level
actn %>%
  select(iteration, moa, starts_with("moa")) %>% 
  distinct() %>%
  # Examine SD across classes for each iteration
  group_by(iteration) %>% 
  summarise_at(vars(starts_with("moa_")), sd) %>%
  # Summarise results across all iterations
  summarise_at(vars(starts_with("moa_")), mean)

saveRDS(actn, file = "scratch_data/sim2_interactn_opts.Rds")
save(rheum_final, file = "scratch_Data/sim2_interaction_opts_ordering.Rds")
