#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=60:00:00
#PBS -l walltime=72:00:00
/usr/bin/Rscript simuln/02b_run_inla_models_drug.R 9 hi atc5_0.25_trial_0.25_drug_0.25 &>> simuln/sim1/hi/output.txt
