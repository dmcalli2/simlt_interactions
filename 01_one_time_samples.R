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
colnames(drug) <- paste("drug", colnames(drug),  sep = "_")
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

actn <- bind_cols(atc5, trial, drug) %>% 
  arrange(atc_5, drug, nct_id) %>% 
  mutate(iteration = rep(1:1000, 161)) %>% 
    select(atc_5, drug, nct_id, iteration, everything())
saveRDS(actn, file = "scratch_data/interactn_opts.Rds")
save(diabetes_final, file = "scratch_Data/interaction_opts_ordering.Rds")
  # inner_join(trial)
# 
# CombineVarn <- function(component = "cept", choice = list(atc5 = "0.25",
#                                                        drug = "0.25",
#                                                        trial = "0.25")){
#    atc5[[component]][diabetes_final_smpls$atc_5,  choice[["atc5"]]] +
#    drug[[component]][diabetes_final_smpls$drug, choice[["drug"]]] +
#   trial[[component]][diabetes_final_smpls$nct_id, choice[["trial"]]]
# }
# diabetes_final_smpls$cept <- CombineVarn()
# diabetes_final_smpls$como <- CombineVarn("como", list(atc5 = "0.5",
#                                                             drug = "0.5",
#                                                             trial = "0.5"))
# diabetes_final_smpls$allc <- CombineVarn("allc", list(atc5 = "0.15",
#                                                             drug = "0.2",
#                                                             trial = "0.25"))
# diabetes_final_smpls$actn <- CombineVarn("actn", list(atc5 = "0.05",
#                                                             drug = "0.1",
#                                                             trial = "0.15"))
# 
# ## Add in wider group effect to difference around this to get each effect
# diabetes_final_smpls <- within(diabetes_final_smpls, {
#   cept = cept + main_scen$cept[1]
#   como = como + main_scen$como[1]
#   allc = allc + main_scen$allc[1]
#   actn = actn + main_scen$actn[1]
# })
# 
# ############## DO NEED TO ADD TOGETHER EFFECTS TO GET OBSERVED GROUP MEANS!!!
# diabetes_final_smpls_obs <- within(diabetes_final_smpls, {
#   cept = cept 
#   como = cept + como
#   allc = cept + allc
#   actn = actn + como + allc - cept
# })
# 
# ### Calcualate error terms based on SD
# # Standardised sd is approx 1 for change in HbA1c
# comorbidity_prev <- 0.2
# sd <- 1
# ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
# ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5
# 
# names(ncomo_se) <- diabetes_final$nct_id
# names(ycomo_se) <- diabetes_final$nct_id
# 
# ncomo_prec <- 1/ncomo_se^2
# ycomo_prec <- 1/ycomo_se^2
# 
# ## Add these into both dataframes
# diabetes_final$ncomo_se <- ncomo_se
# diabetes_final$ycomo_se <- ycomo_se
# 
# diabetes_final_smpls_obs$ncomo_se   <- ncomo_se[diabetes_final_smpls_obs$nct_id]
# diabetes_final_smpls_obs$ycomo_se   <- ycomo_se[diabetes_final_smpls_obs$nct_id]
# diabetes_final_smpls_obs$ncomo_prec <- ncomo_prec[diabetes_final_smpls_obs$nct_id]
# diabetes_final_smpls_obs$ycomo_prec <- ycomo_prec[diabetes_final_smpls_obs$nct_id]
# 
