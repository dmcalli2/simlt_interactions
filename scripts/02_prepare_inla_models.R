# 02_create_interaction_datasets.R
library(INLA)

## Read in diabetes trials
load("Data/metadata_for_simulation.Rdata")

# Create effect estimates interactions with variation at trial, drug and class level
# Read in dataset with different interactions
# Each variable is for a different component
diabetes <- readRDS("scratch_data/interactn_opts.Rds")
diabetes <- as.data.frame(diabetes)

# extract vectors describing trial, drug and atc5 level interactions
atc5s  <- names(diabetes)[substr(names(diabetes), 1, 5) == "atc5_"]
trials <- names(diabetes)[substr(names(diabetes), 1, 6) == "trial_"]
drugs <- names(diabetes)[substr(names(diabetes), 1, 5) == "drug_"]

# Choose fewer options to reduce number need to run
atc5s <- atc5s[c(1,3,5)]
trials <- trials[c(1,3,5)]
drugs <- drugs[c(1,3,5)]

# Create matrix of combinations
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

# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

# Calculate standard error for groups with and without comorbidity
ncomo_se = sd/ ((1-comorbidity_prev) * diabetes_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * diabetes_final$n_per_grp)^0.5

# Calculate SE for interaction, same for all
# We multiple by two because the standard error is the same in the treatment and placebo arms
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


## Select only diabetes variables need for each analysis
diabetes <- diabetes [ , c("atc_5", "drug", "nct_id", "iteration")]
save(my_data, myform_nested2, diabetes, res, file = "data/for_inla.Rdata")

### Create scripts to run on HPCC
for(scenario in res_names) {
  con <- file(description =  paste0("unix_scripts/",scenario, ".sh"), open = "wb")
  top <- c("#!/bin/bash",
          "#PBS -l nodes=1:ppn=1:centos6",
          "#PBS -l cput=2:00:00")
  
  act <- paste("/usr/bin/Rscript simuln/02b_run_inla_models.R",
 ## act <- paste("/usr/bin/Rscript simuln/02c_run_inla_class_level.R",
                            scenario,
          ##      "> /export/home/dma24j/run.output", sep = " ")
                  "&>> simuln/output.txt", sep = " ")
  readr::write_lines(c(top, act), con)
  close(con)
}
# Metascript to run scripts
# CAn run up to 50 at a time on short list
# all take about 30 mins
con <- file(description =  "unix_scripts/metascript.sh", open = "wb")
metascript <- paste0("qsub simuln/", res_names, ".sh")
readr::write_lines(metascript, con)
close(con)



