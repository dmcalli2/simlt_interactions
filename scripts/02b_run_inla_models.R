#02b_run_inla_model
library(INLA)
INLA:::inla.dynload.workaround() 
with_args <- TRUE

#Como_prevs
como_prev <- c("hi")
#como_prev <- c("std")
#como_prev <- c("lo")

load(file = paste0("data/sim1/",como_prev,"/for_inla.Rdata"))
############ From now on putty, pass with args
## Loop through 6 scenarios, this will take approximately 3 hours

argsd <- commandArgs(trailingOnly=TRUE)

print(argsd)

choose_scenario <- ifelse(with_args, argsd[3] , "path_0.25_moa_0.05_trial_0.05_drug_0.05")

print(choose_scenario)

# Add in wider drug group level effect to a chosen variation scenario
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
  saveRDS(scenario, file = paste0("simuln/sim1/",como_prev,"/output/sim1_",argsd[1],"_",argsd[2],"_",choose_scenario, ".rds" ))
