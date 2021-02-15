# 02_create_interaction_datasets.R
library(INLA)
library(magrittr)
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

#Como_prevs
como_prev <- c("hi")
#como_prev <- c("std")
#como_prev <- c("lo")

comorbidity_prev <- ifelse(como_prev == "hi", 0.5, NA) #cardiovascular disease
comorbidity_prev <- ifelse(como_prev == "std", 0.2, comorbidity_prev)
comorbidity_prev <- ifelse(como_prev == "lo", 0.1, comorbidity_prev) #copd/repiratory condition

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

# Make linear combinations
d1 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(1,rep(NA,24-1))) # drug level
d2 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,1),1,rep(NA,24-2))) # drug level
d3 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,2),1,rep(NA,24-3))) # drug level
d4 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,3),1,rep(NA,24-4))) # drug level
d5 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,4),1,rep(NA,24-5))) # drug level
d6 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, 1, NA), mydrug = c(rep(NA,5),1,rep(NA,24-6))) # drug level
d7 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,6),1,rep(NA,24-7))) # drug level
d8 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,7),1,rep(NA,24-8))) # drug level
d9 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA), mydrug = c(rep(NA,8),1,rep(NA,24-9))) # drug level
d10 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, 1, NA, NA, NA, NA, NA), mydrug = c(rep(NA,9),1,rep(NA,24-10))) # drug level
d11 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,10),1,rep(NA,24-11))) # drug level
d12 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,11),1,rep(NA,24-12))) # drug level
d13 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,12),1,rep(NA,24-13))) # drug level
d14 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(1, NA, NA, NA, NA, NA, NA), mydrug = c(rep(NA,13),1,rep(NA,24-14))) # drug level
d15 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,14),1,rep(NA,24-15))) # drug level
d16 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,15),1,rep(NA,24-16))) # drug level
d17 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,16),1,rep(NA,24-17))) # drug level
d18 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, NA, NA, 1), mydrug = c(rep(NA,17),1,rep(NA,24-18))) # drug level
d19 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,18),1,rep(NA,24-19))) # drug level
d20 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, 1, NA, NA, NA, NA), mydrug = c(rep(NA,19),1,rep(NA,24-20))) # drug level
d21 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,20),1,rep(NA,24-21))) # drug level
d22 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,21),1,rep(NA,24-22))) # drug level
d23 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, NA, 1, NA, NA), mydrug = c(rep(NA,22),1,rep(NA,24-23))) # drug level
d24 <- inla.make.lincomb(myatc4 = 1, myatc5 = c(NA, NA, NA, 1, NA, NA, NA), mydrug = c(rep(NA,23),1) )# drug level

## Write drug-only model
myform_nested3 <- y ~ -1 + mydrug + 
  f(trial, model = "iid", 
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
save(my_data, myform_nested2,myform_nested3, diabetes, res,
     d1,d2,d3,d4,d5,d6,
     d7,d8,d9,d10,d11,
     d12,d13,d14,d15,
     d16,d17,d18,d19,
     d20,d21,d22,d23,d24, file = paste0("data/sim1/",como_prev,"/for_inla.Rdata"), version = 2)

diabetes$res <- res[, 'atc5_0.05_trial_0.05_drug_0.05'] + -0.1
my_data$y <-  diabetes$res[diabetes$iteration == 1]

my_data_for_table <- cbind(diabetes_final, my_data) %>%
  dplyr::select(nct_id, atc_5, drug, n_per_grp, y_prec, y) %>%
  dplyr::mutate(sd = sqrt(y_prec^-1) ) 

#write.csv(my_data_for_table, file= "./tables/sim1my_data.csv")  

### Create scripts to run on HPCC
scenarios <- c("atc5_0.05_trial_0.05_drug_0.05",
               "atc5_0.15_trial_0.05_drug_0.05",
               "atc5_0.25_trial_0.05_drug_0.05",
               "atc5_0.05_trial_0.15_drug_0.05",
               "atc5_0.05_trial_0.25_drug_0.05",
               "atc5_0.05_trial_0.05_drug_0.15",
               "atc5_0.05_trial_0.05_drug_0.25",
               "atc5_0.15_trial_0.15_drug_0.15",
               "atc5_0.25_trial_0.25_drug_0.25")

modelvect <- c(drug = "simuln/02b_run_inla_models_drug.R")#, 
               #scen = "simuln/02b_run_inla_models_scenario.R")
count <- 0
for(scenario in scenarios) {
  for(modelchoose in names(modelvect)) {
  count <- count + 1
  con <- file(description =  paste0("unix_scripts/sim1/",como_prev,"/",
                                    count,"_",
                                    como_prev,"_",
                                    scenario,"_", 
                                    modelchoose, ".sh"), open = "wb")
  top <- c("#!/bin/bash",
           "#PBS -l nodes=1:ppn=1:centos6",
           "#PBS -l cput=60:00:00",
       "#PBS -l walltime=72:00:00") ## setting walltime at >4hrs is necessary for low prevalence runs, but should be avoided for others as it bumps up into the long queue
  
  act <- paste("/usr/bin/Rscript", paste0(modelvect[modelchoose]),
               ## act <- paste("/usr/bin/Rscript simuln/2_02c_run_inla_class_level.R",
               count,  como_prev ,scenario,
               ##      "> /export/home/dma24j/run.output", sep = " ")
               "&>>", paste0("simuln/sim1/",como_prev,"/output.txt"), sep = " ")
  readr::write_lines(c(top, act), con)
  close(con)
}}

# Metascript to run scripts
# CAn run up to 50 at a time on short list
# all take about 30 mins
# a <- as.character(c(seq(1,length(res_names),1)))
# res_names2 <- paste(a,como_prev,res_names, sep="_")

a <- list.files(paste0("unix_scripts/sim1/", como_prev))
a <- setdiff(a, "metascript.sh")
# a <- a[stringr::str_detect(a, "drug\\.sh$")]
metascript <- paste0("qsub simuln/sim1/",como_prev,"/", a)
con <- file(description =  paste0("unix_scripts/sim1/",como_prev,"/metascript.sh"), open = "wb")
readr::write_lines(metascript, con)
close(con)

