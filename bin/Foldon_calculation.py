import numpy as np
import Domains
import argparse
from multiprocessing import Pool
import subprocess, random, os

# Modify the following lines according to your NUPACK installation:
nupack_path = os.environ['HOME'] + '/nupack3.2.2/build/bin'
nupack_env = {'NUPACKHOME' : os.environ['HOME'] + '/nupack3.2.2'}
subopt_gap = 0.99

# Change following routines for other environments:
L_init = 5  # Initiation unit
dL = 5  # elongation unit (also means CG unit)
MULTI_PROCESS = 192


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


def nupack_subopt(sequence, T=37, gap=subopt_gap):
    # Use NUPACK to calculate the minimum-free-energy secondary structure of a sequence
    # NOTE: Returns a secondary structure in the (((.))) notation
    seq = sequence
    # print(seq)
    rint = int(random.random() * 1.e9)
    tmp = nupack_path + '/tmp/%d' % rint
    with open(tmp + '.in', 'w') as f:
        f.write("%s\n" % seq)
        f.write("%f\n" % gap)
    subprocess.call([nupack_path + '/subopt', '-T', str(T), tmp])  # , env=nupack_env) Using system path env
    sss = []
    with open(tmp + '.subopt', 'r') as f:
        flag = False
        for line in f:
            if len(line) > 1 and all(c in '(.)' for c in line.strip()):
                ss = line.strip()
                sss.append(ss)
    os.remove(tmp + '.in')
    os.remove(tmp + '.subopt')
    return sss


def save_foldon(l_bound, r_bound, ss, f):
    f.write(f'{l_bound} {r_bound} {ss} \n')
    return True


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence prefix (one line input)")
    parser.add_argument('--path', type=str, default="foldons.dat", help="path to store foldons")
    parser.add_argument('--mfe-only', action='store_true', help="Use mfe structures only")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    #   NOTE: Initiation [create active population]
    all_foldons = Domains.FoldonCollection()
    sequence_length = len(full_sequence)
    sequence_length = 300
    foldons = open(clargs.path, 'w+')
    if clargs.mfe_only:
        for current_length in range(L_init, sequence_length+dL, dL):
            l_bounds = np.arange(0, current_length, dL)
            multi_pool = Pool()
            new_foldons_ss = list(multi_pool.map(nupack_mfe, [full_sequence[l_bound:current_length] for l_bound in l_bounds]))
            multi_pool.close()
            multi_pool.join()
            for i in range(len(l_bounds)):
                save_foldon(l_bounds[i], current_length, new_foldons_ss[i], foldons)
            foldons.flush()
    else:
        for current_length in range(L_init, sequence_length, dL):
            l_bounds = np.arange(0, current_length, dL)
            multi_pool = Pool()
            new_foldons_sss = list(multi_pool.map(nupack_subopt, [full_sequence[l_bound:current_length] for l_bound in l_bounds]))
            multi_pool.close()
            multi_pool.join()
            for i in range(len(l_bounds)):
                for ss in new_foldons_sss[i]:
                    save_foldon(l_bounds[i], current_length, ss, foldons)
            foldons.flush()
        current_length = sequence_length
        l_bounds = np.arange(0, current_length, dL)
        multi_pool = Pool()
        new_foldons_sss = list(multi_pool.map(nupack_subopt, [full_sequence[l_bound:current_length] for l_bound in l_bounds]))
        multi_pool.close()
        multi_pool.join()
        for i in range(len(l_bounds)):
            for ss in new_foldons_sss[i]:
                save_foldon(l_bounds[i], current_length, ss, foldons)
        foldons.flush()

    foldons.close()
    exit()
