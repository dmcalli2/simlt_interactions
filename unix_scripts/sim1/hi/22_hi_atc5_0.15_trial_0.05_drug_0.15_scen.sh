#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=48:00:00
#PBS -l walltime=50:00:00
/usr/bin/Rscript simuln/02b_run_inla_models_scenario.R 22 hi atc5_0.15_trial_0.05_drug_0.15 &>> simuln/sim1/hi/output.txt
