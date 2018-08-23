#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=2:00:00
/usr/bin/Rscript simuln/02b_run_inla_models.R atc5_0.05_trial_0.05_drug_0.15 &>> simuln/output.txt
