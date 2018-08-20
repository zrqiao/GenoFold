from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import argparse, math, random, gzip, pickle, types
from collections import defaultdict
import os
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


def localss_population_processing(input_prefix):
    local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
    if not os.path.exists(input_prefix):
        os.makedirs(input_prefix)
    with open(input_prefix + '.dat', 'r+') as folding_input:
        sss = [(x.split()[0], np.float(x.split()[1]))
               for x in folding_input.readlines()]
        for ss in sss:
            if ss[0].startswith('#'):
                time = ss[1]
            else:
                # print(ss)
                if len(ss[0]) >= SD_end:
                    SD_ss = ss[0][SD_start:SD_end]
                    local_structure_collection_data[SD_ss][time] += ss[1]
    with open(input_prefix + '/local_population' + '.dat', 'w+') as local_output:
        for local_ss in local_structure_collection_data.keys():
            local_output.write(local_ss + '\n')
            local_output.write(' '.join(map(str, local_structure_collection_data[local_ss].keys())) + '\n')
            local_output.write(' '.join(map(str, local_structure_collection_data[local_ss].values())) + '\n')


def data_ploting(ax_punbound, input_prefix, label, start_index, end_index):
    data_punbound = defaultdict(np.float)
    if not os.path.exists(input_prefix):
        os.makedirs(input_prefix)
    if not os.path.exists(input_prefix + '/p_unbound'):
        os.makedirs(input_prefix + '/p_unbound')
    f = open(input_prefix + f'/p_unbound/base{start_index}_{end_index}.dat', 'w')
    with open(input_prefix + '.dat', 'r+') as folding_input:
        sss = [(x.split()[0], np.float(x.split()[1]))
               for x in folding_input.readlines()]
        for ss in sss:
            if ss[0].startswith('#'):
                time = ss[1]
            else:
                # print(ss)
                if len(ss[0]) >= SD_end:
                    target_ss = ss[0][start_index:end_index]
                    data_punbound[time] += ss[1] * target_ss.count('.') / (end_index - start_index)
        data_plot = np.array([list(data_punbound.keys()), list(data_punbound.values())])

    ax_punbound.plot(data_plot[0], data_plot[1], label=label)
    for d in data_punbound.items():
        f.write(f'{d[0]}  {d[1]}\n')
    f.close()


def data_ploting_equ(ax_punbound, input_prefix, label, start_index, end_index):
    data_punbound = defaultdict(np.float)
    if not os.path.exists(input_prefix):
        os.makedirs(input_prefix)
    if not os.path.exists(input_prefix + '/p_unbound'):
        os.makedirs(input_prefix + '/p_unbound')
    f = open(input_prefix + f'/p_unbound/base{start_index}_{end_index}.dat', 'w')
    with open('folA_WT/summary_pairs.dat', 'r+') as folding_input:  # NOTE: Format here need to be unified!
        pairs_data = [list(map(np.float, x.split())) for x in folding_input.readlines()]
        for dat in pairs_data:
            if len(dat) >= end_index:
                time = len(dat)-1  # NOTE: should be len(dat) for later version
                data_punbound[time] += np.sum(dat[start_index:end_index]) / (end_index - start_index)
        data_plot = np.array([list(data_punbound.keys()), list(data_punbound.values())])

    ax_punbound.plot(data_plot[0], data_plot[1], label=label)
    for d in data_punbound.items():
        f.write(f'{d[0]}  {d[1]}\n')
    f.close()


if __name__ == '__main__':

    # plt.style.use('bmh')
    # mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize'] = 20
    mpl.rcParams['axes.titleweight'] = 20
    params = {
        'axes.labelsize': 15,
        'legend.fontsize': 12.5,
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'text.usetex': False,
        'figure.figsize': [15, 9]
    }
    mpl.rcParams.update(params)
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    parser.add_argument('--working-path', type=str, default='.', help="Path to store outputs")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open(clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    PATH = clargs.working_path
    # Start IO
    fig = plt.figure()
    fig.add_axes()
    NUM_COLORS = 6
    cm = plt.get_cmap('Set1')
    ax_punbound = fig.add_subplot(111)
    ax_punbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

    # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {PATH}')
    ax_punbound.set_title(f'{PATH}: '+r'Average SD $p_{unbound}$')
    ax_punbound.set_xlabel('Transcription time')
    ax_punbound.set_ylabel(r'$p_{unbound}$')
    ax_punbound.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    # ax_localpop.set_yscale('log')
    # ax_localpop.set_xscale('log')
    # ax_localpop.set_ylim(1e-5, 1.5)
    # ax_localpop.set_ylim(0.0, 1.1)

    for e_k in range(1, 15, 3):
        k = 1*10**e_k
        print('k= %.2g'%k)
        prefix = PATH + '/k' + '%.2g' % k
        label = r'$k/k_T$ = ' + '%.2g' % k
        # localss_population_processing(prefix)
        data_ploting(ax_punbound, prefix, label, SD_start, SD_end)

    print('k= inf')
    data = defaultdict(np.float)
    local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
    prefix = PATH + '/k' + 'inf'
    label = r'$k/k_T$ = ' + 'inf'
    # localss_population_processing(prefix)
    data_ploting(ax_punbound, prefix, label, SD_start, SD_end)

    prefix = PATH + '/equilibrium'  # Need to be copied here
    label = 'Equilibrium'
    data_ploting_equ(ax_punbound, prefix, label, SD_start, SD_end)  # Another format
    ax_punbound.legend(loc='best')
    # fig.tight_layout()
    plt.show()

    fig.savefig(PATH + f'/p_unbound_SD_k_tuning.eps')

    for base_position in range(0, 7):  # note: Relative to SD_start
        base_gene_position = base_position-14
        fig = plt.figure()
        fig.add_axes()
        ax_punbound = fig.add_subplot(111)
        ax_punbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {clargs.sequence}')
        ax_punbound.set_title(f'{PATH}: base [{base_gene_position}] '+r'$p_{unbound}$:')
        # ax_punbound.set_title(f'{clargs.sequence} base [{base_position}] '+r'$p_{unbound}$:')
        ax_punbound.set_xlabel('Transcription time')
        ax_punbound.set_ylabel(r'$p_{unbound}$')
        # ax_punbound.set_yscale('log')
        # ax_localpop.set_xscale('log')
        # ax_punbound.set_ylim(1e-4, 1.1)
        ax_punbound.grid(axis='y', color="0.9", linestyle='-', linewidth=1)

        for e_k in range(1, 15, 3):
            k = 1 * 10 ** e_k
            print('k= %.2g' % k)
            prefix = PATH + '/k' + '%.2g' % k
            label = r'$k/k_T$ = ' + '%.2g' % k
            data_ploting(ax_punbound, prefix, label, SD_start+base_position, SD_start+base_position+1)

        print('k= inf')
        data = defaultdict(np.float)
        local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
        prefix = PATH + '/k' + 'inf'
        label = r'$k/k_T$ = ' + 'inf'
        data_ploting(ax_punbound, prefix, label, SD_start+base_position, SD_start+base_position+1)

        prefix = PATH + '/equilibrium'
        label = 'Equilibrium'
        data_ploting_equ(ax_punbound, prefix, label, SD_start+base_position, SD_start+base_position+1)

        ax_punbound.legend(loc='best')
        # fig.tight_layout()
        plt.show()

        fig.savefig(PATH + f'/p_unbound_base[{base_gene_position}]_k_tuning.eps')
exit()
