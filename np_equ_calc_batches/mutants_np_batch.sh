#!/bin/bash
#SBATCH -n 192                # Number of cores
#SBATCH -N 12                # Ensure that all cores are on one machine
#SBATCH -t 15-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem=240000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid

module load python
cat mutants_list.in | while read prefix
do
    echo $prefix
    seq=$(head -n 1 sequences/$prefix.in)
    echo ${#seq}
    mkdir $prefix/segments
    python np_scripts/segments_align.py $prefix
    sbatch np_scripts/get_segment_pairs.sh $prefix
done


# monitor output (need formatted string)

