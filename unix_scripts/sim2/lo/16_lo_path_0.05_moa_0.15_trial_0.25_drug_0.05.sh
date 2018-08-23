#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=4:00:00
#PBS -l walltime=8:00:00
/usr/bin/Rscript simuln/sim2/lo/2_02b_run_inla_models.R 16 lo path_0.05_moa_0.15_trial_0.25_drug_0.05 &>> simuln/sim2/lo/output.txt
