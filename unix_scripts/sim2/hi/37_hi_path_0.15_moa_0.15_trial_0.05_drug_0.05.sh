#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=4:00:00
#PBS -l walltime=8:00:00
/usr/bin/Rscript simuln/sim2/hi/2_02b_run_inla_models.R 37 hi path_0.15_moa_0.15_trial_0.05_drug_0.05 &>> simuln/sim2/hi/output.txt
