#02b_run_inla_model
library(INLA)
# library(tidyverse)
INLA:::inla.dynload.workaround() 
with_args <- TRUE
############ From now on putty, pass with args
argsd <- commandArgs(trailingOnly=TRUE)
print(argsd)

como_prev <- ifelse(with_args, argsd[2] ,"std")

load(file = paste0("data/sim1/",como_prev,"/for_inla.Rdata"))

## Loop through 6 scenarios, this will take approximately 3 hours

choose_scenario <- ifelse(with_args, argsd[3] , "atc5_0.15_trial_0.15_drug_0.05")

print(choose_scenario)
print(como_prev)

# Add in wider drug group level effect to a chosen variation scenario
 diabetes$res <- res[, choose_scenario] + -0.1

# Loop through each iteration
 scenario <- list()
 drugs <- list()
 for (iter in 1:1000){
# browser()
   print(iter)
    ## Add values for specific iteration
    my_data$y <-  diabetes$res[diabetes$iteration == iter]

    ## Run model, trial within drug within ATC5 class within ATC4 class
    mod1_nested2 <- inla(myform_nested2,
                 data = my_data,
                 # Add linear combinations to estimate drug-class
                 lincomb = c(d1= d1,d2= d2, d3= d3, d4= d4, d5= d5, d6= d6,
                             d7= d7,d8= d8,d9= d9,d10= d10, d11= d11,
                             d12= d12, d13= d13, d14= d14, d15= d15,
                             d16= d16, d17= d17, d18= d18, d19= d19,
                             d20= d20, d21= d21, d22= d22, d23= d23, d24= d24),
                 # Likelihood distribution
                 family = "gaussian",
                 # Fix likelihood hyperpars as the data is fixed with known precision.
                 control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                 # Likelihood precisions
                 scale = my_data$y_prec,
                 # Prior distribution for "fixed" effects - really for mu_mu
                 control.fixed = list(mean = 0, prec = 0.25),
                 # Optionally compute DIC
                 control.compute = list(config=TRUE,dic = TRUE, waic = TRUE),
                 verbose = FALSE)
    # Storing iterations in a single list
    scenario[[iter]] <- summary(mod1_nested2)

 }
 
  dir.create(paste0("simuln/sim1/",como_prev,"/output/"), showWarnings = FALSE)
  # Saving each list as a data file
  saveRDS(scenario, file = paste0("simuln/sim1/",
                                  como_prev,
                                  "/output/sim1_",argsd[1],"_",argsd[2],"_",
                                  choose_scenario, "_full.rds" ))
  

