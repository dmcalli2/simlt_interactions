# 02_create_interaction_datasets.R
# Create effect estimates interactions with variation at trial, drug and class level
diabetes <- readRDS("scratch_data/interactn_opts.Rds")

count <- 0
res <- matrix(nrow = nrow(diabetes), ncol = 5^3)
res_names <- vector(length = 5^3)

for (i in c("atc5_0.05", "atc5_0.1", "atc5_0.15", "atc5_0.2", "atc5_0.25")){
  for(j in c("trial_0.05", "trial_0.1", "trial_0.15", "trial_0.2", "trial_0.25")){
    for(k in c("drug_0.05", "drug_0.1", "drug_0.15", "drug_0.2", "drug_0.25")){
      count <- count + 1
    res[, count] <- 
       diabetes[ , i] +
       diabetes[ , j] +
       diabetes[ , k] 
    res_names[count] <- paste(i, j, k, sep = "_")
      }
  }
}
colnames(res) <- res_names
rownames(res) <- paste(diabetes$atc_5, diabetes$drug, diabetes$nct_id, diabetes$iteration,
                       sep = "_")

# Add in se term for interaction
comorbidity_prev <- 0.2
sd <- 1
# calculate SE for comorbidity adn non-comorbidity group (same for placebo and treatment)
## Read in diabetes trials
load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")

# Remove single alpha-glucosidase inhibitor trial
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5
# Calculate SE for interaction, same for all
inter_prec <- 1/(2*ncomo_se^2 + 2*ycomo_se^2)

library(INLA)

## Write model
myform_nested2 <- y ~ -1 + myatc4 + 
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01)))) +
  f(myatc5, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.01))))

## Make part of model matrix which is identical for all iterations
my_drug_n <- as.numeric(as.factor(diabetes_final$drug))
study_id_n <- as.numeric(as.factor(as.character(diabetes_final$nct_id)))
atc5_n <- as.numeric(as.factor(as.character(diabetes_final$atc_5)))

## Create dataset which is consistent for all iterations
my_data <- data.frame(y_prec = inter_prec, 
                      trial = study_id_n,
                      myatc4 = 1,
                      myatc5 = atc5_n,
                      mydrug = my_drug_n)

## Loop through 6 scenarios, this will take approximately 3 hours
scenarios <- c("atc5_0.05_trial_0.05_drug_0.05",
                  "atc5_0.1_trial_0.1_drug_0.1",
                  "atc5_0.25_trial_0.25_drug_0.25",
                  "atc5_0.1_trial_0.1_drug_0.25",
                  "atc5_0.1_trial_0.25_drug_0.1",
                  "atc5_0.25_trial_0.1_drug_0.1")

for (choose_scenario in scenarios) {
print(choose_scenario)
# Add in main effect to a chosen variation scenario
 diabetes$res <- res[, choose_scenario] + -0.2

# Loop through each iteration
 scenario <- vector(length = 250, mode = "list")
 for (iter in 1:250){
      
    ## Add values for specific iteration
    my_data$y <-  diabetes$res[diabetes$iteration == iter]
    
    ## Run model, trial within drug within ATC5 class within ATC4 class
    mod1_nested2 <- inla(myform_nested2, 
                 data = my_data,
                 # Add linear combinations to estimate drug-class
                 # lincomb = dc_all,
                 # Likelihood distribution
                 family = "gaussian",
                 # Fix likelihood hyperpars as the data is fixed with known precision.
                 control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                 # Likelihood precisions
                 scale = my_data$y_prec,
                 # Prior distribution for "fixed" effects - really for mu_mu
                 control.fixed = list(mean = 0, prec = 0.25),
                 # Optionally compute DIC
                 verbose = FALSE)
    # Storing iterations in a single list
    scenario[[iter]] <- summary(mod1_nested2)
 }
 # Saving each list as a data file
 saveRDS(scenario, file = paste0("simulate_interactions/scenario_",choose_scenario, ".rds" ))
}
# 
# 
# fxd <- lapply(scenario, function(x) x$fixed)
# fxd <- do.call(rbind, fxd)
