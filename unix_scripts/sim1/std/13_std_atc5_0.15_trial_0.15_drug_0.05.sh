#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=4:00:00
#PBS -l walltime=6:00:00
/usr/bin/Rscript simuln/sim1/std/02b_run_inla_models.R 13 std atc5_0.15_trial_0.15_drug_0.05 &>> simuln/sim1/std/output.txt
