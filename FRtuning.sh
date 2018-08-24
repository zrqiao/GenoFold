#!/bin/bash
#SBATCH -n 48                # Number of cores
#SBATCH -N 3                # Ensure that all cores are on one machine
#SBATCH -t 0-48:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem-per-cpu=16000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
sequence=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACATATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAA

module load python
mut=$1
for ((i=0;i<45;i=i+1))
do
    echo -e "Simulating with k = $i"
    python bin/GenoFold.py --working-path ./$mut/CG5 --k 1e$i --foldons-path $mut/foldons_ext_CG5.dat --pool-size 25 --CG-length 5 sequences/$mut &
    sleep 1m
done
# python bin/GenoFold.py --working-path ./$mut/CG5 --stationary --foldons-path $mut/foldons_ext_CG5.dat --pool-size 25 --CG-length 5 sequences/$mut
sleep 500m
# monitor output (need formatted string)

