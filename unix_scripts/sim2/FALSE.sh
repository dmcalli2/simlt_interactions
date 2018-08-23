#!/bin/bash
#PBS -l nodes=1:ppn=1:centos6
#PBS -l cput=2:00:00
/usr/bin/Rscript simuln/02b_run_inla_models.R FALSE > /export/home/dma24j/run.output
