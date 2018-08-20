#!/bin/bash
#SBATCH -n 16                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-96:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakgpu   # Partition to submit to
#SBATCH --mem=4000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
sequence=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACATATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAA

module load python

#for ((i=1;i<=14;i=i+1))
#do
#    echo -e "Simulating with k = $i"
#    python bin/GenoFold.py --working-path ./folA_WT/ext_fds_CG5_pool50 --k 1e$i --foldons-path folA_WT/foldons_subopt_CG.dat --pool-size 50 --CG-length 5 folA_WT &
#    sleep 1m
#done
python bin/GenoFold.py --working-path ./folA_WT/ext_fds_CG5_pool1000 --stationary --foldons-path folA_WT/foldons_subopt_CG.dat --pool-size 1000 --CG-length 5 folA_WT

# monitor output (need formatted string)

