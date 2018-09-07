#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-96:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakgpu   # Partition to submit to
#SBATCH --mem=16000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o truncated_pairs.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e truncated_pairs.err  # File to which STDERR will be written, %j inserts jobid

mut=$1
export NUPACKHOME=/n/home02/jarrenqiao/nupack3.2.2
module load python
seq=$(head -n 1 sequences/$mut.in)
for ((i=1;i<=${#seq};i=i+1))
do
    echo -e "Analyzing segment $i"
    # ../build/bin/pfunc -T 37 segments/sample_$i > segments/sample_$i.pfunc
    $NUPACKHOME/build/bin/mfe -T 37 $mut/segments/sample_$i
    # ../build/bin/prob -T 37 segments/sample_$i > segments/sample_$i.prob
    $NUPACKHOME/build/bin/pairs -T 37 $mut/segments/sample_$i
done
python np_scripts/segments_result_parsing.py $mut
