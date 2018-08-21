#!/bin/bash
#SBATCH -n 64                # Number of cores
#SBATCH -N 4                # Ensure that all cores are on one machine
#SBATCH -t 0-128:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem=160000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
5UTR=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACAT
module load python

mkdir $prefix
python bin/Foldon_calculation.py sequences/$prefix --path $prefix/foldons_ext_CG5.dat


# monitor output (need formatted string)

