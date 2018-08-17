#!/bin/bash
#SBATCH -n 32                # Number of cores
#SBATCH -N 2                # Ensure that all cores are on one machine
#SBATCH -t 0-96:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem=160000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
sequence=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACATATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAA

module load python

for ((i=1;i<=15;i=i+1))
do
    echo -e "Simulating with k = $i"
    python GenoFold.py --working-path ./folA_WT/ext_fds --k 1e$i --foldons-path folA_WT/foldons.dat folA_WT &
    sleep 1m
done
python GenoFold.py --working-path ./folA_WT/ext_fds --stationary --foldons-path folA_WT/foldons.dat folA_WT

# monitor output (need formatted string)

