# 02_create_interaction_datasets.R
library(INLA)

## Read in inflamm/GI trials
load("Data/rheum_final.Rdata")

# Create effect estimates interactions with variation at trial, drug and class level
# Read in dataset with different interactions
# Each variable is for a different component
rheum <- readRDS("scratch_data/sim2_interactn_opts.Rds")
rheum <- as.data.frame(rheum)

## Select only one MoA - Interleukin-17A Antagonists & Receptor Interactions 
ilk17a <- rheum %>% 
  filter(moa == "Interleukin-17A Antagonists & Receptor Interactions")

# extract vectors describing trial and drug level interactions only
trials <- names(ilk17a)[substr(names(ilk17a), 1, 6) == "trial_"]
drugs <- names(ilk17a)[substr(names(ilk17a), 1, 5) == "drug_"]

# Choose fewer options to reduce number need to run
# atc5s <- atc5s[c(1,3,5)]
trials <- trials[c(1,3,5)]
drugs <- drugs[c(1,3,5)]

# Create matrix of combinations
count <- 0
res <- matrix(nrow = nrow(ilk17a), 
                          ncol = length(trials) * length(drugs))
res_names <- vector(length =  length(trials) * length(drugs))

for(j in trials){
  for(k in drugs){
    count <- count + 1
  res[, count] <- 
     ilk17a[ , j] +
     ilk17a[ , k]
  res_names[count] <- paste(j, k, sep = "_")
    }
}

colnames(res) <- res_names
rownames(res) <- paste(ilk17a$drug, ilk17a$nct_id, ilk17a$iteration,
                       sep = "_")

# Add in se term for interaction

#Como_prevs
como_prev <- c("std")
#como_prev <- c("std")
#como_prev <- c("lo")

comorbidity_prev <- ifelse(como_prev == "hi", 0.5, NA) #cardiovascular disease
comorbidity_prev <- ifelse(como_prev == "std", 0.2, comorbidity_prev)
comorbidity_prev <- ifelse(como_prev == "lo", 0.1, comorbidity_prev) #copd/repiratory condition

sd <- 1
# calculate SE for comorbidity adn non-comorbidity group (same for placebo and treatment)

# Select only one MOA
ilk17a_final <-  
  subset(rheum_final, moa == "Interleukin-17A Antagonists & Receptor Interactions")

ilk17a_final <- ilk17a_final[with(ilk17a_final, order(drug, nct_id)),]

# Calculate standard error for groups with and without comorbidity
ncomo_se = sd/ ((1-comorbidity_prev) * ilk17a_final$n_per_grp)^0.5
ycomo_se = sd/ (   comorbidity_prev  * ilk17a_final$n_per_grp)^0.5

# Calculate SE for interaction, same for all
# We multiple by two because the standard error is the same in the treatment and placebo arms
inter_prec <- 1/(2*ncomo_se^2 + 2*ycomo_se^2)

## Write model
myform_nested2 <- y ~ -1 + wdg +    ## wdg = wider (/widest) drug grouping (top of hierarchy)
  f(trial, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.1)))) +
  f(mydrug, model = "iid", 
    hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 0.1)))) #+
  # f(myatc5, model = "iid", 
    # hyper = list(prec = list(prior = "logtnormal", param = c(mean = 0, prec = 1))))

## Make part of model matrix which is identical for all iterations
my_drug_n <- as.numeric(as.factor(ilk17a_final$drug))
study_id_n <- as.numeric(as.factor(as.character(ilk17a_final$nct_id)))

## Create dataset which is consistent for all iterations
my_data <- data.frame(y_prec = inter_prec, 
                      trial = study_id_n,
                      wdg = 1,
                      mydrug = my_drug_n)


## Select only diabetes variables need for each analysis
ilk17a <- ilk17a [ , c("drug", "nct_id", "iteration")]
save(my_data, myform_nested2, ilk17a, res, file = paste0("data/sim2/",como_prev,"/one_class_for_inla.Rdata"))

### Create scripts to run on HPCC
count <- 0
for(scenario in res_names) {
  count <- count + 1
  con <- file(description =  paste0("unix_scripts/sim2/",como_prev,"/oneclass_",count,"_",como_prev,"_",scenario, ".sh"), open = "wb")
  top <- c("#!/bin/bash",
           "#PBS -l nodes=1:ppn=1:centos6",
           "#PBS -l cput=2:00:00",
       "#PBS -l walltime=4:00:00") ## setting walltime at >4hrs is necessary for low prevalence runs, but should be avoided for others as it bumps up into the long queue
  
  act <- paste("/usr/bin/Rscript", paste0("simuln/sim2/",como_prev,"/02b_run_inla_models_one_class.R"),
               ## act <- paste("/usr/bin/Rscript simuln/2_02c_run_inla_class_level.R",
               count,  como_prev ,scenario,
               ##      "> /export/home/dma24j/run.output", sep = " ")
               "&>>", paste0("simuln/sim2/",como_prev,"/output_one_class.txt"), sep = " ")
  readr::write_lines(c(top, act), con)
  close(con)
}

# Metascript to run scripts
# CAn run up to 50 at a time on short list
# all take about 30 mins
a <- as.character(c(seq(1,length(res_names),1)))
res_names2 <- paste(a,como_prev,res_names, sep="_")

con <- file(description =  paste0("unix_scripts/sim2/",como_prev,"/oneclass_metascript.sh"), open = "wb")
metascript <- paste0("qsub simuln/sim2/",como_prev,"/oneclass_", res_names2, ".sh")
readr::write_lines(metascript, con)
close(con)




