#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=2:00:00
#PBS -l walltime=4:00:00
/usr/bin/Rscript simuln/sim1/std/02b_run_inla_models_one_class.R 5 std trial_0.15_drug_0.15 &>> simuln/sim1/std/output_one_class.txt
