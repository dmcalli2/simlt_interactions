#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=4:00:00
#PBS -l walltime=6:00:00
/usr/bin/Rscript simuln/sim1/lo/02b_run_inla_models.R 3 lo atc5_0.05_trial_0.05_drug_0.25 &>> simuln/sim1/lo/output.txt
