#02b_run_inla_model
library(INLA)
INLA:::inla.dynload.workaround() 
inputty <- TRUE

############ From now on putty, pass with args
## Loop through 6 scenarios, this will take approximately 3 hours
scenarios <- c("atc5_0.05_trial_0.05_drug_0.05",
                  "atc5_0.1_trial_0.1_drug_0.1",
                  "atc5_0.25_trial_0.25_drug_0.25",
                  "atc5_0.1_trial_0.1_drug_0.25",
                  "atc5_0.1_trial_0.25_drug_0.1",
                  "atc5_0.25_trial_0.1_drug_0.1")
## Loop through 1 scenario, this will take approximately 3 hours
scenarios <- c("atc5_0.05_trial_0.05_drug_0.05")#,
                  # "atc5_0.1_trial_0.1_drug_0.1",
                  # "atc5_0.25_trial_0.25_drug_0.25",
                  # "atc5_0.1_trial_0.1_drug_0.25",
                  # "atc5_0.1_trial_0.25_drug_0.1",
                  # "atc5_0.25_trial_0.1_drug_0.1")
scenarios

for (choose_scenario in scenarios) {
print(choose_scenario)
# Add in main effect to a chosen variation scenario
 diabetes$res <- res[, choose_scenario] + -0.1

# Loop through each iteration
 scenario <- vector(length = 250, mode = "list")
 for (iter in 1:1000){
      
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