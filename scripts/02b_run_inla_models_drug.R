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

choose_scenario <- ifelse(with_args, argsd[3] , "atc5_0.05_trial_0.05_drug_0.05")

print(choose_scenario)
print(como_prev)

# Add in wider drug group level effect to a chosen variation scenario
 diabetes$res <- res[, choose_scenario] + -0.1

# Loop through each iteration
 scenario <- list()
 drugs <- list()
 for (iter in 1:1000){

    ## Add values for specific iteration
    my_data$y <-  diabetes$res[diabetes$iteration == iter]
    # trial_count <- tapply(my_data$mydrug, my_data$mydrug, length)
    # trial_count <- which(trial_count == 1)
    # my_data <- my_data[!my_data$mydrug %in% trial_count, ]

    drug_models <- list()
    
    for(j in 1:length(unique(my_data$mydrug)))
    {
      #browser()
      # Reduce down to include data only from individual drug
      sel_drug_mydata <- my_data[my_data$mydrug==j, ]
      sel_drug_mydata$mydrug <- 1 #Have to recode to 1 in order to be consistent across drugs
      
      ## Run model, trial within drug 
      mod1_nested3 <- inla(myform_nested3, 
                           data = sel_drug_mydata,
                           #data = my_data_test,            #Use instead if testing 
                           # Add linear combinations to estimate drug-class
                           #lincomb = d1,
                           #lincom = do.call(c, dc_tests),  #Use instead if testing 
                           # Likelihood distribution
                           family = "gaussian",
                           # Fix likelihood hyperpars as the data is fixed with known precision.
                           control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                           # Likelihood precisions
                           # scale = my_data_test$y_prec,   #Use instead if testing 
                           scale = sel_drug_mydata$y_prec,
                           # Prior distribution for "fixed" effects - really for mu_mu
                           control.fixed = list(mean = 0, prec = 0.25),
                           # Optionally compute DIC
                           verbose = FALSE,
                           control.compute = list(config=TRUE,dic = TRUE, waic = TRUE),
                           control.inla = list())
      
      drug_models[[j]] <- summary(mod1_nested3)
    }
    drugs[[iter]] <- drug_models
 }
  dir.create(paste0("simuln/sim1/",como_prev,"/output/"), showWarnings = FALSE)
  # Saving each list as a data file
  saveRDS(drugs, file = paste0("simuln/sim1/"
                               ,como_prev,"/output/sim1_",
                               argsd[1],"_",argsd[2],"_",
                               choose_scenario, "_drug.rds" ))
