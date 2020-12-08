#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=48:00:00
#PBS -l walltime=50:00:00
/usr/bin/Rscript simuln/02b_run_inla_models_scenario.R 50 lo atc5_0.25_trial_0.25_drug_0.05 &>> simuln/sim1/lo/output.txt
