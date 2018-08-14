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
SD_start, SD_end = 21, 28
equi_p_unbound = [0.0414220, 0.0612670, 0.0839040, 0.9764600, 0.9300200, 0.0861740, 0.2976000]


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def local_plot(ax_localpop, local_input_path, label):
    with open(local_input_path + '.dat', 'r+') as local_input:
        # ax_localpop.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {clargs.sequence}')
        ax_localpop.set_title(f'SD Local folding population $k/k_T$ = {label}')
        ax_localpop.set_xlabel('Transcription time')
        ax_localpop.set_ylabel('Population fraction')
        ax_localpop.set_xlim(0, 515)
        print(f'SD Local folding population $k/k_T$ = {label}')
        data_raw = local_input.readlines()
        ss_num = int(len(data_raw) / 3)
        for i in range(ss_num):
            ss = data_raw[i * 3].rstrip('\n')
            times = list(map(np.float, data_raw[i * 3 + 1].split()))
            populations = list(map(np.float, data_raw[i * 3 + 2].split()))
            ax_localpop.plot(times, populations, label=ss)
    ax_localpop.legend(loc='best')


if __name__ == '__main__':

    plt.style.use('ggplot')
    fig = plt.figure(figsize=(10, 10))
    # colors = [plt.cm.jet(lt) for lt in range(0, 8)]
    fig.add_axes()

    # mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize'] = 10
    mpl.rcParams['axes.titleweight'] = 10
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    fig = plt.figure(figsize=(12, 12))
    fig.add_axes()

    for e_k in range(1, 16, 1):
        ax_localpop = fig.add_subplot(4, 4, e_k)
        k = 1 * 10 ** e_k
        local_input_path = clargs.sequence + '_local_population_k' + '%e' % k
        label = '%.2g' % k
        local_plot(ax_localpop, local_input_path, label)

    ax_localpop = fig.add_subplot(4, 4, 16)

    local_input_path = clargs.sequence + '_local_population_k' + 'inf'
    label = 'inf'
    local_plot(ax_localpop, local_input_path, label)

    fig.tight_layout()
    plt.show()

    fig.savefig(clargs.sequence + '_local_population_evolution_summary.eps')

    exit()
