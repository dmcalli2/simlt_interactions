# 02_create_interaction_datasets.R
library(INLA)

## Read in inflamm/GI trials
load("Data/rheum_final.Rdata")

# Create effect estimates interactions with variation at trial, drug and class level
# Read in dataset with different interactions
# Each variable is for a different component
rheum <- readRDS("scratch_data/sim2_interactn_opts.Rds")
rheum <- as.data.frame(rheum)

# extract vectors describing trial, drug and atc5 level interactions

paths  <- names(rheum)[substr(names(rheum), 1, 5) == "path_"]
moas  <- names(rheum)[substr(names(rheum), 1, 4) == "moa_"]
trials <- names(rheum)[substr(names(rheum), 1, 6) == "trial_"]
drugs <- names(rheum)[substr(names(rheum), 1, 5) == "drug_"]

# Choose fewer options to reduce number need to run
paths <- paths[c(3)]    # Reduce number of scenarios by sampling only one level of variation at highest level
moas <- moas[c(1,3,5)]
trials <- trials[c(1,3,5)]  
drugs <- drugs[c(1,3,5)]

# Create matrix of combinations
count <- 0
res <- matrix(nrow = nrow(rheum), ncol = length(paths) *length(moas) * length(trials) * length(drugs))
res_names <- vector(length = length(paths) *length(moas) * length(trials) * length(drugs))

for (h in paths){
  for(i in moas){
    for(j in trials){
      for(k in drugs){
      count <- count + 1
    res[, count] <- 
       rheum[ , h] + 
       rheum[ , i] +
       rheum[ , j] +
       rheum[ , k]
    res_names[count] <- paste(h, i, j, k, sep = "_")
      }
    }  
  }
}
colnames(res) <- res_names
rownames(res) <- paste(rheum$brd_drug_pth, rheum$moa, rheum$drug, rheum$nct_id, rheum$iteration,
                       sep = "_")

# Add in se term for interaction
comorbidity_prev <- 0.2

load("data/outcome_smrs_for_simulation.Rdata")
print(c(das_smrs_sd, ibdq_smrs_sd))

# Values of 0.14 and 0.17 from 'standardized' versions of DAS & IBDQ
# If correct, I can't see value of distinguishing - will set to 0.2 here
# but could presumably use 1 as diabetes (not sure how HB1c was scaled) - check
# this with DM

sd <- 0.2       ### This should come from real studies and depend on indication (i.e., on outcome IBDQ, DAS)
      
# calculate SE for comorbidity adn non-comorbidity group (same for placebo and treatment)

rheum_final <- rheum_final[with(rheum_final, order(atc_5, drug, nct_id)),]

# Calculate standard error for groups with and without comorbidity
ncomo_se = sd/ ((1-comorbidity_prev) * rheum_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * rheum_final$n_per_grp)^0.5

# Calculate SE for interaction, same for all
# We multiple by two because the standard error is the same in the treatment and placebo arms
inter_prec <- 1/(2*ncomo_se^2 + 2*ycomo_se^2)

## Write model
myform_nested2 <- y ~ -1 + wdg +    ## wdg = wider (/widest) drug grouping (top of hierarchy)
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mymoa, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1)))) +
  f(mypath, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1))))

## Make part of model matrix which is identical for all iterations
my_drug_n <- as.numeric(as.factor(rheum_final$drug))
study_id_n <- as.numeric(as.factor(as.character(rheum_final$nct_id)))
moa_n <- as.numeric(as.factor(as.character(rheum_final$moa)))
path_n <- as.numeric(as.factor(as.character(rheum_final$brd_drug_pth)))

## Create dataset which is consistent for all iterations
my_data <- data.frame(y_prec = inter_prec, 
                      trial = study_id_n,
                      wdg = 1,            ## wdg = wider (/widest) drug grouping (top of hierarchy)
                      mymoa = moa_n,
                      mydrug = my_drug_n,
                      mypath = path_n)


## Select only rheum variables need for each analysis
rheum <- rheum [ , c("brd_drug_pth", "moa", "drug", "nct_id", "iteration")]
save(my_data, myform_nested2, rheum, res, file = "data/sim2_for_inla.Rdata")

### Create scripts to run on HPCC
for(scenario in res_names) {
  con <- file(description =  paste0("unix_scripts/sim2/",scenario, ".sh"), open = "wb")
  top <- c("#!/bin/bash",
          "#PBS -l nodes=1:ppn=1:centos6",
          "#PBS -l cput=2:00:00")
  
  act <- paste("/usr/bin/Rscript simuln/02b_run_inla_models.R",
               ## act <- paste("/usr/bin/Rscript simuln/02c_run_inla_class_level.R",
               scenario,
               ##      "> /export/home/dma24j/run.output", sep = " ")
               "&>> simuln/output_sim2.txt", sep = " ")
  readr::write_lines(c(top, act), con)
  close(con)
}
# Metascript to run scripts
# CAn run up to 50 at a time on short list
# all take about 30 mins
con <- file(description =  "unix_scripts/sim2/metascript.sh", open = "wb")
metascript <- paste0("qsub simuln/", res_names, ".sh")
readr::write_lines(metascript, con)
close(con)



