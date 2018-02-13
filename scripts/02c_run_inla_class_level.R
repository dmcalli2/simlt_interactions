#2c run_INLA model class_level analysis

#02b_run_inla_model
library(INLA)
# INLA:::inla.dynload.workaround() 
with_args <- FALSE
load(file = "data/for_inla.Rdata")
############ From now on putty, pass with args
## Loop through 6 scenarios, this will take approximately 3 hours

choose_scenario <- "atc5_0.05_trial_0.05_drug_0.05"
if(with_args) choose_scenario <- commandArgs(trailingOnly=TRUE)

# Add in main effect to a chosen variation scenario
 diabetes$res <- res[, choose_scenario] + -0.1
 diabetes$res[diabetes$atc_5 %in% c("A10BA","A10BB")] <- NA

 
 # Mkae linear combination
 dc1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(1, NA, NA, NA, NA, NA, NA, NA)) 
 dc2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA, NA)) 
 dc_all <- c(dc1,dc2)
 names(dc_all) <- paste0("dc", 1:2)
# Loop through each iteration
 # scenario <- vector(length = 1, mode = "list")
 # for (iter in 1:1){
      iter <- 1
    ## Add values for specific iteration
    my_data$y <-  diabetes$res[diabetes$iteration == iter]
    
    ## Run model, trial within drug within ATC5 class within ATC4 class
    mod1_nested2 <- inla(myform_nested2, 
                 data = my_data,
                 # Add linear combinations to estimate drug-class
                 lincomb = dc_all,
                 # Likelihood distribution
                 family = "gaussian",
                 # Fix likelihood hyperpars as the data is fixed with known precision.
                 control.family = list(hyper = list(prec = list(fixed = TRUE, initial = 0))),
                 # Likelihood precisions
                 scale = my_data$y_prec,
                 # Prior distribution for "fixed" effects - really for mu_mu
                 control.fixed = list(mean = 0, prec = 0.25),
                 # Optionally compute DIC
                 verbose = FALSE,
                 control.compute = list(config=TRUE))
    # Storing iterations in a single list
    # scenario[[iter]] <- summary(mod1_nested2)
     summary(mod1_nested2)
     a <- mod1_nested2$summary.random
     mod1_nested2$summary.fixed
     a$myatc5
     summary(mod1_nested2)
     mod1_nested2$summary.lincomb.derived
    
     
     smpls <- inla.posterior.sample(100000, mod1_nested2)
     smpls <- lapply(smpls, function(x) x$latent[c("myatc5:1", "myatc5:2", "myatc4"),])
     smpls <- do.call(rbind, smpls)
     smpls <- smpls[,1:2] + smpls[,c(3,3)]

     t(apply(smpls, 2, function(x) {
        c(mean = mean(x),
             median = median(x),
             lci = quantile(x, probs = c(0.025)),
             uci = quantile(x, probs = c(0.975)),
             sd = sd(x))}))
      
     plot(density(smpls[,1]), col = "blue", xlim = c(-2, 2))
     mydnorm <- function(x) dnorm(x, mean = mean(smpls[,1]), sd = sd(smpls[,1]))
     curve(mydnorm, from = -2, to = 2, add = TRUE, col = "black")
     
     mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
     fitted_t <- MASS::fitdistr(smpls[,1], mydt, list(m = -0.1, s = 0.5, df = 3),
                                lower = -5)
     res <- fitted_t$estimate
     mydt_fitted <- function(x) mydt(x, res["m"], res["s"], res["df"])
     curve(mydt_fitted, from = -2, to = 2, add = TRUE, col = "red")
     
    ## T distribution fits data well (as per BUGS)

    

    
 # }
 # Saving each list as a data file
 # saveRDS(scenario, file = paste0("scenario_class_lvl",choose_scenario, ".rds" ))