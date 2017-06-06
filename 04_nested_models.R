## 04_nested_models.R

drugs_select <- "Antidiabetic"
source('02_arrange_data_arrays.R') 

library(rjags)
#load.module('glm')
library(tidyverse)
library(stringr)

#### Model with study nested inside drug-class ----
modelstringnested_class <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], prec_matrix[, , Dstudy])
        for(i in 1:Ncoef) { # 8 independent priors for first 8 coefficients
        mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        mu[9, Dstudy]  ~ dnorm(dep_class[Ddrug],  mu9_prec[Ddrug])
        mu[10, Dstudy] ~ dnorm(pain_class[Ddrug], mu10_prec[Ddrug])
      } # end study
      # Priors for each condition
      }# end condition
 # Priors for each drug class
      dep_class[Ddrug]  ~ dnorm(dep, 0.0001)           # mean dep:alloc
      pain_class[Ddrug] ~ dnorm(pain, 0.0001)           # mean pain:alloc
      mu9_sd[Ddrug]  ~ dnorm(0,0.0001)T(0,)      # sd dep:alloc
      mu10_sd[Ddrug] ~ dnorm(0,0.00001)T(0,)      # sd pain:alloc
      mu9_prec[Ddrug] <- 1/mu9_sd[Ddrug]^2   # prec dep:alloc
      mu10_prec[Ddrug] <- 1/mu10_sd[Ddrug]^2 # prec pain:alloc
  }# end drug class
  dep  ~ dnorm(0, 0.0001)           # mean dep:alloc
  pain ~ dnorm(0, 0.0001)           # mean pain:alloc
}# model
"
writeLines(modelstringnested_class, con= "jags/modelstringnested_class.txt")

#### Model with study nested inside drug-class, informative ----
modelstringnested_class_inform <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], prec_matrix[, , Dstudy])
        for(i in 1:Ncoef) { # 8 independent priors for first 8 coefficients
        mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        mu[9, Dstudy]  ~ dnorm(dep_class[Ddrug],  mu9_prec[Ddrug])
        mu[10, Dstudy] ~ dnorm(pain_class[Ddrug], mu10_prec[Ddrug])
      } # end study
      # Priors for each condition
      }# end condition
 # Priors for each drug class
      dep_class[Ddrug]  ~ dnorm(dep, 0.25)           # mean dep:alloc
      pain_class[Ddrug] ~ dnorm(pain, 0.25)           # mean pain:alloc
      mu9_sd[Ddrug]  ~ dnorm(0,0.00001)T(0,)      # sd dep:alloc
      mu10_sd[Ddrug] ~ dnorm(0,0.00001)T(0,)      # sd pain:alloc
      mu9_prec[Ddrug] <- 1/mu9_sd[Ddrug]^2   # prec dep:alloc
      mu10_prec[Ddrug] <- 1/mu10_sd[Ddrug]^2 # prec pain:alloc
  }# end drug class
  dep  ~ dnorm(0, 0.25)           # mean dep:alloc
  pain ~ dnorm(0, 0.25)           # mean pain:alloc
}# model
"
writeLines(modelstringnested_class_inform, con= "jags/modelstringnested_class_inform.txt")


#### Model with study only, pooled ----
modelstring_pooled <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], prec_matrix[, , Dstudy])
        for(i in 1:Ncoef) { # 8 independent priors for first 8 coefficients
        mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        mu[9, Dstudy]  ~ dnorm(dep,  mu9_prec)
        mu[10, Dstudy] ~ dnorm(pain, mu10_prec)
      } # end study
      # Priors for each condition
      }# end condition
 # Priors for each drug class
  }# end drug class
      dep  ~ dnorm(0, 0.0001)           # mean dep:alloc
      pain ~ dnorm(0, 0.0001)           # mean pain:alloc
      mu9_sd  ~ dnorm(0,0.00001)T(0,)      # sd dep:alloc
      mu10_sd ~ dnorm(0,0.00001)T(0,)      # sd pain:alloc
      mu9_prec <- 1/mu9_sd^2   # prec dep:alloc
      mu10_prec <- 1/mu10_sd^2 # prec pain:alloc
}# model
"
writeLines(modelstring_pooled, con= "jags/modelstring_pooled.txt")

#### Model with study only, fixed ----
modelstring_fixed <- "
model{
  for(Ddrug in 1:Ndrug) {
    for(Dcondition in Sdrug[Ddrug]:(Sdrug[Ddrug+1]-1)){
      for(Dstudy in Scondition[Dcondition]:(Scondition[Dcondition+1]-1)){
        coef[, Dstudy] ~ dmnorm(mu[,Dstudy], prec_matrix[, , Dstudy])
        for(i in 1:Ncoef) { # 8 independent priors for first 8 coefficients
        mu[i, Dstudy] ~ dnorm(0, 0.0001)
        }
        mu[9, Dstudy]  <- dep
        mu[10, Dstudy] <- pain
      } # end study
      # Priors for each condition
      }# end condition
 # Priors for each drug class
  }# end drug class
      dep  ~ dnorm(0, 0.0001)           # mean dep:alloc
      pain ~ dnorm(0, 0.0001)           # mean pain:alloc
}# model
"
writeLines(modelstring_fixed, con= "jags/modelstring_fixed.txt")


#### Run models ----
# Turn list of coefficients into a matrix of coefficients
# Where each column is a study

## Loop through effect estimates

map(seq_along(res_rag), function (effect_number) {
  
  ## Select effect estimate
  myrag <- res_rag[[effect_number]]
  
  ## Chose outcome type, or loop through
  outcomes <- list(con = myrag$con, cat = myrag$cat) # con continuous and cat categorical
  
  ## Chose modeltype, or loop through
  modeltypes <- c(fixed = "jags/modelstring_fixed.txt", 
                 pooled = "jags/modelstring_pooled.txt",
                 dc_nest = "jags/modelstringnested_class.txt",
                 dc_nest_inform = "jags/modelstringnested_class_inform.txt")
  
  params_monitor <- list(fixed = c('dep', 'pain'),
                         pooled = c('dep', 'pain'),
                         dc_nest = c('dep_class', 'pain_class', 'dep', 'pain', 'mu9_sd', 'mu10_sd'),
                         dc_nest_inform = c('dep_class', 'pain_class', 'dep', 'pain', 'mu9_sd', 'mu10_sd'))
  
  for(outcome in names(outcomes)){
    # Extract coefficients into matrix, each row is a trial
    coef <-  map(outcomes[[outcome]], ~ .x$coef)
    coef <- sapply(coef, identity, simplify="matrix")
    print(dim(coef))
    
    # Turn list of precision matrices into array of precision matrices
    # Where each z is a study
    prec_matrix <-  map(outcomes[[outcome]], ~ .x$prec_matrix) 
    prec_matrix <- sapply(prec_matrix, identity, simplify="array")
    prec_matrix <- round(prec_matrix, 3)
    print(dim(prec_matrix))
  
    for(modeltype in names(modeltypes)){
        # Run Jags Model    
        jags <- jags.model(modeltypes[modeltype],
                           data = list (Ndrug = myrag$Ndrug,
                                        Sdrug = myrag$Sdrug,
                                        Scondition = myrag$Scondition,
                                        coef = coef,
                                        prec_matrix = prec_matrix,
                                        Ncoef = 8 # note 8 of 10 coefficients are independent for each study
                           ),
                           n.chains = 2,
                           n.adapt = 1000)
        diag <- coda.samples(jags, c("mu", params_monitor[[modeltype]]), 10000)
        update(jags, 10000) # Burn in
        
        data(LINE)
        LINE$recompile()
        LINE.out <- coda.samples(jags,
                                 params_monitor[[modeltype]],
                                 20000)
        print(paste0(outcome, "_", modeltype, ".Rds"))
        saveRDS(summary(LINE.out), file = paste0("model_summaries/", outcome, "_", modeltype, "_effect", effect_number, ".Rds"))
        saveRDS(diag, file = paste0("mcmc_chains/", outcome, "_", modeltype, "_effect", effect_number, ".Rds"))
    } ## End of models loop
  } ## End of outcomes loop
})
