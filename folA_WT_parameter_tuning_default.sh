#!/bin/bash
#SBATCH -n 16                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-96:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shakhnovich   # Partition to submit to
#SBATCH --mem=40000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
sequence=ATACCCGTTTTTTGGGCTAACAGGAGGAATTACATATGATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAA

module load python

for ((i=1;i<=14;i=i+3))
do
    echo -e "Simulating with k = $i"
    python bin/GenoFold.py --working-path ./folA_WT/ext_fds_pool100 --k 1e$i --foldons-path folA_WT/foldons_ext.dat --pool-size 100 --CG-length 1 folA_WT &

    # ../build/bin/mfe -T 37 -dangles some segments/sample_$i
    sleep 1m
done
python bin/GenoFold.py --working-path ./folA_WT/ext_fds_pool100 --stationary --foldons-path folA_WT/foldons_ext.dat --pool-size 100 --CG-length 1 folA_WT


# monitor output (need formatted string)

