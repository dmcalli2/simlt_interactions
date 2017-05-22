# 00 Read in data
library(readr)
library(stringr)
library(dplyr)
library(tidyr)

## Read in csv with selected trials
trials <- read_csv(file = "../export_r/Summary conditions and treatments reviewed_Main.csv", col_types = paste(rep("c",17),collapse = "") )
trials$drug_group_allocation [trials$drug_group_allocation == "Haemostatsis"] <- "Haemostasis"

## Select non discontinued drug types - 171 trials
trials <- filter (trials, ! bnf_group %in% c("Discontinued", "0"))

## Read in number of participants
clinical_study <- read_delim("../AACT201409Annual_Pipe_delimited_txt/clinical_study_noclob_nolf.txt", delim = "|")  
names(clinical_study) <- tolower(names(clinical_study))

clinical_study <- filter (clinical_study, nct_id %in% trials$nct_id)

## Sample of NCTID by sponsor
set.seed(1234)
sponsor_sample <- group_by(clinical_study, source)
sponsor_sample <- sample_n(sponsor_sample, size = 1)
# Most had data on websites (75%)

## Join
trials <- inner_join(trials, clinical_study[ , c("nct_id", "enrollment")])
rm(clinical_study)
gc()

save (trials, file = "Dataset_of_trials.Rdata")
