#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 4-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem=12000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o batch_mutants_pop_CG5_pool50.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e batch_mutants_pop_CG5_pool50.err  # File to which STDERR will be written, %j inserts jobid
sequence=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACATATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAA

module load python
mut=$1
i=$2
echo -e "Simulating with k = $i"
python bin/GenoFold.py --working-path ./$mut/CG5_pool40 --k 1e$i --foldons-path $mut/foldons_ext_CG5.dat --pool-size 40 --CG-length 5 sequences/$mut
# python bin/GenoFold.py --working-path ./$mut/CG5 --stationary --foldons-path $mut/foldons_ext_CG5.dat --pool-size 25 --CG-length 5 sequences/$mut
# monitor output (need formatted string)

