#02b_run_inla_model
library(INLA)
INLA:::inla.dynload.workaround() 
with_args <- TRUE
load(file = "data/sim2withpath_for_inla.Rdata")
############ From now on putty, pass with args
## Loop through 6 scenarios, this will take approximately 3 hours

choose_scenario <- "path_0.15_moa_0.05_trial_0.05_drug_0.05"
if(with_args) choose_scenario <- commandArgs(trailingOnly=TRUE)

print(choose_scenario)

# Add in wider drug group level effect to a chosen variation scenario
rheum$res <- res[, choose_scenario] + -0.1

# Loop through each iteration
 scenario <- vector(length = 250, mode = "list")
 for (iter in 1:1000){

    ## Add values for specific iteration
    my_data$y <-  rheum$res[rheum$iteration == iter]
    
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
  saveRDS(scenario, file = paste0("sim2_scenario_",choose_scenario, ".rds" ))
