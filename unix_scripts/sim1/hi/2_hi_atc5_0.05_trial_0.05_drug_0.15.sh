#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=4:00:00
#PBS -l walltime=6:00:00
/usr/bin/Rscript simuln/sim1/hi/02b_run_inla_models.R 2 hi atc5_0.05_trial_0.05_drug_0.15 &>> simuln/sim1/hi/output.txt
