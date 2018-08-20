from difflib import SequenceMatcher
import numpy as np
from bin import Domains, nupack_functions
from multiprocessing import Pool
import argparse
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


def pfunc_map(dat):
    time = dat[0]
    seq = dat[1]
    print(time)
    pfunc = nupack_functions.nupack_pfunc(seq, 37)
    return (time, pfunc)


def boltzmann_map(dat):
    time = dat[0]
    seq = dat[1]
    ss = dat[2]
    # print(ss)
    G = nupack_functions.nupack_ss_free_energy(seq, ss, 37)
    bolz = Domains.rate(G, 1)
    return (time, ss, bolz)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    print('k = inf')
    # data = defaultdict(np.float)
    local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
    pfuncs = defaultdict(np.float)
    preprocess_data1 = []  # time, seq, ss
    preprocess_data2 = []  # time, seq
    # f = open(clargs.sequence + '_p_unbound_inf_nupack' + '.dat', 'w')
    with open(clargs.sequence + '_stationary' + '.dat', 'r+') as folding_input:

        sss = [(x.split()[0], np.float(x.split()[1]))
               for x in folding_input.readlines()]

        for ss in sss:
            if ss[0].startswith('#'):
                print(ss)
                time = ss[1]
                seq = full_sequence[:int(5*np.ceil(time/5))]
                preprocess_data2.append([time, seq])
            else:
                if len(ss[0]) >= SD_end:
                    preprocess_data1.append([time, seq, ss[0]])
                    # norm_c[time] += nupack_prob
                    # data[time] += ss[1] * SD_ss.count('.') / (SD_end - SD_start)

        # data_p = np.array([list(data.keys()), list(data.values())])

    pool2 = Pool()
    data2 = pool2.map(pfunc_map, preprocess_data2)
    pool2.close()
    pool = Pool()
    data1 = pool.map(boltzmann_map, preprocess_data1)
    pool.close()

    for dat in data2:
        pfuncs[dat[0]] = dat[1]
    for dat in data1:
        time = dat[0]
        ss = dat[1]
        bolz = dat[2]
        nupack_prob = bolz / pfuncs[time]
        SD_ss = ss[SD_start:SD_end]
        local_structure_collection_data[SD_ss][time] += nupack_prob

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
