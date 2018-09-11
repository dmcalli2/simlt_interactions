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
atc5s <- atc5s[c(1,5)]
trials <- trials[c(1,5)]
drugs <- drugs[c(1,5)]

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

# Select only highest and lowest variation scenarios
res <- res[,c('atc5_0.05_trial_0.05_drug_0.05','atc5_0.25_trial_0.25_drug_0.25')]

res_names <- colnames(res)

# Add in se term for interaction

#Como_prevs
como_prev <- c("range")
#como_prev <- c("std")
#como_prev <- c("lo")

if(como_prev == "range") comorbidity_prev <-seq(0.02,0.2,0.02)


sd <- 1
# calculate SE for comorbidity adn non-comorbidity group (same for placebo and treatment)
for(cprev in comorbidity_prev) {
# Remove single alpha-glucosidase inhibitor trial, already out of diabetes file
diabetes_final <-  
  subset(diabetes_final, diabetes_final$atc_5 != "A10BF")
diabetes_final <- diabetes_final[with(diabetes_final, order(atc_5, drug, nct_id)),]

# Calculate standard error for groups with and without comorbidity
ncomo_se = sd/ ((1-cprev) * diabetes_final$n_per_grp)^0.5
ycomo_se = sd/ (   cprev  * diabetes_final$n_per_grp)^0.5

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
save(my_data, myform_nested2, diabetes, res, file = paste0("data/sim1/",como_prev,"/for_inla_prev",cprev,".Rdata"))
}

###Modified to here

### Create scripts to run on HPCC

count <- 0
for(scenario in res_names) {
  for(cprev in comorbidity_prev){
  count <- count + 1
  con <- file(description =  paste0("unix_scripts/sim1/",como_prev,"/",count,"_prev",cprev,"_",scenario, ".sh"), open = "wb")
  top <- c("#!/bin/bash",
           "#PBS -l nodes=1:ppn=1:centos6",
           "#PBS -l cput=4:00:00",
       "#PBS -l walltime=6:00:00") ## setting walltime at >4hrs is necessary for low prevalence runs, but should be avoided for others as it bumps up into the long queue
  
  act <- paste("/usr/bin/Rscript", paste0("simuln/sim1/",como_prev,"/02b_run_inla_models.R"),
               ## act <- paste("/usr/bin/Rscript simuln/2_02c_run_inla_class_level.R",
               count,  cprev ,scenario,
               ##      "> /export/home/dma24j/run.output", sep = " ")
               "&>>", paste0("simuln/sim1/",como_prev,"/output.txt"), sep = " ")
  readr::write_lines(c(top, act), con)
  close(con)
  }
}

# Metascript to run scripts
# CAn run up to 50 at a time on short list
# all take about 30 mins
a <- as.character(c(seq(1,20,1)))
b <- as.character(c(comorbidity_prev,comorbidity_prev))
c <- c(rep(res_names[1],10), rep(res_names[2],10))
res_names2 <- paste0(a,"_prev",b,"_",c)

con <- file(description =  paste0("unix_scripts/sim1/",como_prev,"/metascript.sh"), open = "wb")
metascript <- paste0("qsub simuln/sim1/",como_prev,"/", res_names2, ".sh")
readr::write_lines(metascript, con)
close(con)




