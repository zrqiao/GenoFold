from difflib import SequenceMatcher
import numpy as np
import Domains
import nupack_functions
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import argparse, math, random, gzip, pickle, types
from collections import defaultdict

# Change following routines for other environments:
L_init = 10  # Initiation unit
dL = 10  # elongation unit (also means CG unit)
transcription_time = 0.1
dt = transcription_time * dL  # Folding time for each elongation step (0.1 s/nt)
population_size_limit = 100  # maximum type of strands in the pool
MULTI_PROCESS = 32
SD_start, SD_end = 21, 28
prob_cutoff = 0.001
equi_p_unbound = [0.0414220, 0.0612670, 0.0839040, 0.9764600, 0.9300200, 0.0861740, 0.2976000]


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    print('k= inf')
    # data = defaultdict(np.float)
    local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
    norm_c = defaultdict(np.float)
    # f = open(clargs.sequence + '_p_unbound_inf_nupack' + '.dat', 'w')
    with open(clargs.sequence + '_stationary' + '.dat', 'r+') as folding_input:

        sss = [(x.split()[0], np.float(x.split()[1]))
               for x in folding_input.readlines()]

        for ss in sss:
            if ss[0].startswith('#'):
                print(ss)
                time = ss[1]
                pfunc=0
            else:
                # print(ss)
                if len(ss[0]) >= SD_end:
                    seq = full_sequence[:len(ss[0])]
                    SD_ss = ss[0][SD_start:SD_end]
                    G = nupack_functions.nupack_ss_free_energy(seq, ss[0], 37)
                    if pfunc == 0:
                        pfunc = nupack_functions.nupack_pfunc(seq, 37)
                    nupack_prob = Domains.rate(G, 1)/pfunc
                    # norm_c[time] += nupack_prob
                    # data[time] += ss[1] * SD_ss.count('.') / (SD_end - SD_start)
                    local_structure_collection_data[SD_ss][time] += nupack_prob
        # data_p = np.array([list(data.keys()), list(data.values())])

    with open(clargs.sequence + '_local_population_k' + 'inf_nupack' + '.dat', 'w+') as local_output:
        for local_ss in local_structure_collection_data.keys():
            local_output.write(local_ss + '\n')
            local_output.write(' '.join(map(str, local_structure_collection_data[local_ss].keys())) + '\n')
            local_output.write(' '.join(map(str, local_structure_collection_data[local_ss].values())) + '\n')

    '''
    for d in data.items():
        f.write(f'{d[0]}  {d[1]}\n')
    f.close() 
    '''

exit()
