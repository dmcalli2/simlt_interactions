# 02_create_interaction_datasets.R
library(INLA)
inputty <- FALSE

## Read in diabetes trials
if(inputty == TRUE){
  setwd("simuln")
  load("data/data_for_simulation.Rdata")
} else {
  load("../Trial_identify/clinical_trials_august_2017/scratch_data/data_for_simulation.Rdata")

}

# Create effect estimates interactions with variation at trial, drug and class level
# Read in dataset with different interactions
# Each variable is for a different component
diabetes <- readRDS("scratch_data/interactn_opts.Rds")

# extract vectors describing trial, drug and atc5 level interactions
atc5s  <- names(diabetes)[substr(names(diabetes), 1, 5) == "atc5_"]
trials <- names(diabetes)[substr(names(diabetes), 1, 6) == "trial_"]
drugs <- names(diabetes)[substr(names(diabetes), 1, 5) == "drug_"]

count <- 0
res <- matrix(nrow = nrow(diabetes), ncol = length(atc5s) * length(trials) * length(drugs))
res_names <- vector(length = length(atc5s) * length(trials) * length(drugs))

for (i in atc5s){
  for(j in trials){
    for(k in drugs){
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

# Remove single alpha-glucosidase inhibitor trial
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5
# Calculate SE for interaction, same for all
inter_prec <- 1/(2*ncomo_se^2 + 2*ycomo_se^2)

## Write model
myform_nested2 <- y ~ -1 + myatc4 + 
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(myatc5, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1))))

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


# 

