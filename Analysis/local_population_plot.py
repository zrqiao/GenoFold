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
ddt = 1
dt = transcription_time * dL  # Folding time for each elongation step (0.1 s/nt)
population_size_limit = 100  # maximum type of strands in the pool
MULTI_PROCESS = 32
SD_start, SD_end = 21, 28
km_start = 9
km_end = 37
k_pre = 1e11
km_interval = 1
equi_p_unbound = [0.0414220, 0.0612670, 0.0839040, 0.9764600, 0.9300200, 0.0861740, 0.2976000]


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def local_plot(ax_localpop, local_input_path, label):
    with open(local_input_path + '.dat', 'r+') as local_input:
        # ax_localpop.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {clargs.sequence}')
        ax_localpop.set_title(f'$k_T$ = {label}')
        ax_localpop.set_xlabel('Transcript length')
        ax_localpop.set_ylabel('Population fraction')
        ax_localpop.set_xlim(28, 200)
        # print(f'SD Local folding population k_T = {label}')
        data_raw = local_input.readlines()
        ss_num = int(len(data_raw) / 3)
        sss=[]
        prev_pop = np.zeros(int(520/ddt))
        time_array = np.arange(0, 520, ddt)
        for i in range(ss_num):
            ss = data_raw[i * 3].rstrip('\n')
            times_raw = list(map(np.float, data_raw[i * 3 + 1].split()))
            populations_raw = list(map(np.float, data_raw[i * 3 + 2].split()))
            populations = np.zeros(int(520/ddt))
            for j in range(len(times_raw)):
                populations[int(times_raw[j]/ddt)] = populations_raw[j]

            ax_localpop.bar(time_array, populations, bottom=prev_pop, label=ss, width=ddt)
            prev_pop += populations
    # ax_localpop.legend(loc='best')


if __name__ == '__main__':

    # plt.style.use('ggplot')
    fig = plt.figure(figsize=(27, 20))
    # colors = [plt.cm.jet(lt) for lt in range(0, 8)]
    fig.add_axes()
    # mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize'] = 10
    mpl.rcParams['axes.titleweight'] = 10

    NUM_COLORS = 25
    cm = plt.get_cmap('rainbow_r')

    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    parser.add_argument('--working-path', type=str, default='.', help="Path to store outputs")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    PATH = clargs.working_path
    with open('sequences/'+clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    for e_k in range(km_start, km_end, km_interval):
        ax_localpop = fig.add_subplot(5, 6, int((e_k-km_start)/km_interval)+1)
        ax_localpop.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        k = 1 * 10 ** e_k
        local_input_path = PATH + '/k' + '%.2g' % k + '/local_population'
        label = '%.2g' % (k_pre/k) + ' nt/s'
        local_plot(ax_localpop, local_input_path, label)
    '''
    ax_localpop = fig.add_subplot(2, 3, 6)
    ax_localpop.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

    local_input_path = PATH + '/kinf' + '/local_population'
    label = 'inf'
    local_plot(ax_localpop, local_input_path, label)
    '''
    fig.tight_layout()
    # plt.show()

    fig.savefig(PATH + '/local_population_evolution_summary_bin1.png')
    exit()
