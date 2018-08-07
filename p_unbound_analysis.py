from difflib import SequenceMatcher
import numpy as np
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
SD_start, SD_end = 22, 29


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


if __name__ == '__main__':

    plt.style.use('ggplot')
    fig = plt.figure(figsize=(10, 10))

    # colors = [plt.cm.jet(lt) for lt in range(0, 8)]
    fig.add_axes()

    # mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize'] = 20
    mpl.rcParams['axes.titleweight'] = 15
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    #   NOTE: Initiation [create active population]

    # NOTE: Initiate transcription

    # print('Population size: '+str(active_species_pool.size))

    # Start IO
    fig = plt.figure(figsize=(10, 7.5))
    fig.add_axes()
    gs = gridspec.GridSpec(1, 1, height_ratios=[1])

    #ax_energy = fig.add_subplot(gs[0, 0])
    #ax_energy.set_title('Free Energy')
    # ax_energy.set_xlabel('Subsequence Length', fontsize=12.5)
    #ax_energy.set_ylabel('Free Energy', fontsize=12.5)

    ax_pbound = fig.add_subplot(gs[0, 0])
    ax_pbound.set_title(f'Average p_unbound for SD sequence {clargs.sequence}')
    ax_pbound.set_xlabel('Transcription time', fontsize=12.5)
    ax_pbound.set_ylabel('p_unbound', fontsize=12.5)

    for e_k in range(5, 19):
        k = 1*10**e_k
        print(f'k= {k}')
        data = defaultdict(np.float)
        with open(clargs.sequence + '_k' + '%e' % k + '.dat', 'r') as folding_input:
            f = open(clargs.sequence + '_p_unbound_%e' % k + '.dat', 'w')
            sss = [(x.split()[0].rstrip('\n'), np.float(x.split()[1]))
                   for x in folding_input.readlines() if not x.startswith('#')]

            for ss in sss:
                time = len(ss) * transcription_time
                # print('N=' + str(N) + ' nt')
                # f.write('#N=' + str(N) + ' nt')
                # sim = np.zeros(len(structure))
                SD_ss = ss[SD_start:SD_end]
                data[time] += similar(SD_ss, '.......')
                # print(data[time])
                # for i in range(len(first_sequences)):
                #     temp_cmp = []
                #     for cmpseq in first_sequences:
                #         temp_cmp.append(similar(cmpseq, first_sequences[i]))
                #     sim[i + N] = np.average(np.array(temp_cmp))
                # print(sim)
            data_p = np.array([list(data.keys()), list(data.values())])
            data_p.transpose()

            ax_pbound.plot(data_p, ls='solid', label='k=' + '%.2g' % k + '$s-1')
            for d in data.items():
                f.write(f'{d[0]}  {d[1]}')
            f.close()
    ax_pbound.legend(loc='best')
    fig.tight_layout()
    plt.show()

    fig.savefig('p_unbound_test_k_tuning.eps')