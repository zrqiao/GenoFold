import shlex, subprocess
from difflib import SequenceMatcher
import numpy as np
import Domains
from collections import defaultdict
import nupack_functions
import argparse, math, random, gzip, pickle, types
from multiprocessing import Pool
import copy
import time
import subprocess, random, os, sys

# Modify the following lines according to your NUPACK installation:
nupack_path = os.environ['HOME'] + '/nupack3.2.2/build/bin'
nupack_env = {'NUPACKHOME' : os.environ['HOME'] + '/nupack3.2.2'}

# Change following routines for other environments:
L_init = 1  # Initiation unit
dL = 1  # elongation unit (also means CG unit)
MULTI_PROCESS = 32



def nupack_mfe(sequence, T=37):
    # Use NUPACK to calculate the minimum-free-energy secondary structure of a sequence
    # NOTE: Returns a secondary structure in the (((.))) notation
    seq = sequence
    rint = int(random.random() * 1.e9)
    tmp = nupack_path + '/tmp/%d' % rint
    with open(tmp + '.in', 'w') as f:
        f.write("%s\n" % seq)
    subprocess.call([nupack_path + '/mfe', '-T', str(T), tmp]) #, env=nupack_env) #Using system path env
    with open(tmp + '.mfe', 'r') as f:
        flag = False
        for line in f:
            if len(line) > 1 and all(c in '(.)' for c in line.strip()):
                ss = line.strip()
    os.remove(tmp + '.in')
    os.remove(tmp + '.mfe')
    return ss


def save_foldon(l_bound, r_bound, ss, f):
    f.write(f'{l_bound} {r_bound} {ss} \n')
    return True


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    parser.add_argument('--path', type=str, default="foldons.dat", help="path to store foldons")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    #   NOTE: Initiation [create active population]
    all_foldons = Domains.FoldonCollection()
    sequence_length = len(full_sequence)

    foldons = open(clargs.path, 'w+')

    for current_length in range(L_init, sequence_length+dL, dL):
        l_bounds = np.arange(0, current_length, dL)
        multi_pool = Pool()

        new_foldons_ss = list(multi_pool.map(nupack_mfe, [full_sequence[l_bound:current_length] for l_bound in l_bounds]))

        multi_pool.close()
        multi_pool.join()
        # print(new_foldons_ss)
        for i in range(len(l_bounds)):
            save_foldon(l_bounds[i], current_length, new_foldons_ss[i], foldons)
        foldons.flush()

    foldons.close()
    exit()
