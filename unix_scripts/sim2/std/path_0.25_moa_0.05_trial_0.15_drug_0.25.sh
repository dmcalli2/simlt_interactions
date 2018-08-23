#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=2:00:00
/usr/bin/Rscript simuln/2_02b_run_inla_models.R path_0.25_moa_0.05_trial_0.15_drug_0.25 &>> simuln/output_sim2withpath.txt
