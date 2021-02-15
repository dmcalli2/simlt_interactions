#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=60:00:00
#PBS -l walltime=72:00:00
/usr/bin/Rscript simuln/02b_run_inla_models_drug.R 2 std atc5_0.15_trial_0.05_drug_0.05 &>> simuln/sim1/std/output.txt
