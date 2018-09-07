from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib as mpl
import argparse, math, random, gzip, pickle, types
from collections import defaultdict
import os
import adjustText
# Change following routines for other environments:
L_init = 5  # Initiation unit
dL = 5  # elongation unit (also means CG unit)
transcription_time = 1
dt = transcription_time * dL  # Folding time for each elongation step (0.1 s/nt)
population_size_limit = 100  # maximum type of strands in the pool
MULTI_PROCESS = 32
km_start = 4
km_end = 19
km_interval = 1
SD_start, SD_end = 21, 28
k_pre = 1e11
equi_p_unbound = [0.0414220, 0.0612670, 0.0839040, 0.9764600, 0.9300200, 0.0861740, 0.2976000]
NUM_COLORS = 17

# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


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
                time = ss[1] - 35
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


def data_ploting(ax_punbound, input_prefix, label, start_index, end_index, color_rank):
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
                time = ss[1] - 35
            else:
                # print(ss)
                if len(ss[0]) >= end_index:
                    target_ss = ss[0][start_index:end_index]
                    data_punbound[time] += ss[1] / (end_index - start_index) * target_ss.count('.')
        data_plot = np.array([list(data_punbound.keys()), list(data_punbound.values())])

    ax_punbound.plot(data_plot[0][:-1], data_plot[1][:-1], color=tableau20[color_rank%20])
    tt = plt.text(data_plot[0][-1], data_plot[1][-1], label, fontsize=14, color=tableau20[color_rank%20])
    for d in data_punbound.items():
        f.write(f'{d[0]}  {d[1]}\n')
    f.close()
    return tt


def data_ploting_equ(ax_punbound, mut, input_prefix, label, start_index, end_index, color_rank):
    data_punbound = defaultdict(np.float)
    if not os.path.exists(input_prefix):
        os.makedirs(input_prefix)
    if not os.path.exists(input_prefix + '/p_unbound'):
        os.makedirs(input_prefix + '/p_unbound')
    f = open(input_prefix + f'/p_unbound/base{start_index}_{end_index}.dat', 'w')
    with open(mut + '/summary_pairs.dat', 'r+') as folding_input:  # NOTE: Format here need to be unified!
        pairs_data = [list(map(np.float, x.split())) for x in folding_input.readlines()]
        for dat in pairs_data:
            if len(dat) >= end_index:
                time = len(dat)-35  # NOTE: should be len(dat) for later version
                data_punbound[time] += np.sum(dat[start_index:end_index]) / (end_index - start_index)
        data_plot = np.array([list(data_punbound.keys()), list(data_punbound.values())])

    ax_punbound.plot(data_plot[0][:-1], data_plot[1][:-1], color=tableau20[color_rank%20])
    tt = plt.text(data_plot[0][-1], data_plot[1][-1], label, fontsize=14, color=tableau20[color_rank%20])
    for d in data_punbound.items():
        f.write(f'{d[0]}  {d[1]}\n')
    f.close()
    return tt


if __name__ == '__main__':

    # plt.style.use('bmh')
    # mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize'] = 17
    mpl.rcParams['axes.titleweight'] = 10
    params = {
        'axes.labelsize': 20,
        'legend.fontsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'text.usetex': False,
        'figure.figsize': [12, 9]
    }
    mpl.rcParams.update(params)
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str, help="RNA sequence (one line)")
    parser.add_argument('--working-path', type=str, default='.', help="Path to store outputs")
    # parser.add_argument('--k', type=np.float, default=1., \
    #                     help="pre exponential factor")
    clargs = parser.parse_args()
    with open('sequences/' + clargs.sequence + '.in', 'r') as sequence_file:
        full_sequence = sequence_file.readline().rstrip('\n')

    PATH = clargs.working_path
    mut = clargs.sequence
    # Start IO
    fig = plt.figure()
    fig.add_axes()

    cm = plt.get_cmap('rainbow')
    ax_punbound = fig.add_subplot(111)
    ax_punbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

    # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {PATH}')
    plt.text(250, 1.08, f'{clargs.sequence}: '+r'Average SD $p_{unbound}$', fontsize=17, ha="center")
    ax_punbound.spines["top"].set_visible(False)
    ax_punbound.spines["bottom"].set_visible(False)
    ax_punbound.spines["right"].set_visible(False)
    ax_punbound.spines["left"].set_visible(False)
    ax_punbound.get_xaxis().tick_bottom()
    ax_punbound.get_yaxis().tick_left()
    plt.text(250, -0.08, 'Transcript length', fontsize=14, ha="center", color="0.3")
    ax_punbound.set_ylabel(r'$p_{unbound}$', color="0.3")
    ax_punbound.grid(axis='y', color="0.9", linestyle='--', linewidth=1)
    # ax_localpop.set_yscale('log')
    # ax_localpop.set_xscale('log')
    ax_punbound.set_xlim(-10, 500)
    ax_punbound.set_ylim(0, 1.0)
    color_rank = 0
    labels = []
    for e_k in range(km_start, km_end, km_interval):
        k = 1*10**e_k
        # print('k= %.2g'%k)
        prefix = PATH + '/k' + '%.2g' % k
        label = r'$k_T$ = ' + '%.2g' % (k_pre/k) + r' nt/s'
        try:
            localss_population_processing(prefix)
            labels.append(data_ploting(ax_punbound, prefix, label, SD_start, SD_end, color_rank))
            color_rank += 1
        except:
            continue

    '''
    print('k= inf')
    data = defaultdict(np.float)
    local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
        prefix = PATH + '/k' + 'inf'
    label = r'$k_T/k_f$ = ' + '0' +  r' nt/s'
    localss_population_processing(prefix)
    color_rank += 1
    labels.append(data_ploting(ax_punbound, prefix, label, SD_start, SD_end, color_rank))
    '''

    prefix = PATH + '/equilibrium'  # Need to be copied here
    label = 'Equilibrium Control'
    labels.append(data_ploting_equ(ax_punbound, mut, prefix, label, SD_start, SD_end, color_rank+1))  # Another format

    # ax_punbound.legend(loc='best')
    adjustText.adjust_text(labels, arrowprops=dict(arrowstyle='-', color='0.7'))
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                    labelbottom="on", left="off", right="off", labelleft="on")
    # fig.tight_layout()
    fig.savefig(PATH + f'/p_unbound_SD_k_tuning.png',  bbox_inches="tight")
    # plt.show()

    for base_position in range(0, SD_end-SD_start):  # note: Relative to SD_start
        base_gene_position = base_position+SD_start-35
        fig = plt.figure()
        fig.add_axes()

        ax_punbound = fig.add_subplot(111)
        ax_punbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

        # ax_localpop.set_title(f'Average p_unbound for base[-9](G) {PATH}')
        plt.text(250, 1.08, f'{PATH}: base [{base_gene_position}] '+r'$p_{unbound}$', fontsize=17, ha="center")
        ax_punbound.spines["top"].set_visible(False)
        ax_punbound.spines["bottom"].set_visible(False)
        ax_punbound.spines["right"].set_visible(False)
        ax_punbound.spines["left"].set_visible(False)
        ax_punbound.get_xaxis().tick_bottom()
        ax_punbound.get_yaxis().tick_left()
        plt.text(250, -0.08, 'Transcript length', fontsize=14, ha="center", color="0.3")
        ax_punbound.set_ylabel(r'$p_{unbound}$', color="0.3")
        ax_punbound.grid(axis='y', color="0.9", linestyle='--', linewidth=1)
        # ax_localpop.set_yscale('log')
        # ax_localpop.set_xscale('log')
        ax_punbound.set_xlim(20, 500)
        ax_punbound.set_ylim(0, 1.0)
        color_rank = 0
        labels = []
        for e_k in range(km_start, km_end, km_interval):
            k = 1 * 10 ** e_k
            # print('k= %.2g' % k)
            prefix = PATH + '/k' + '%.2g' % k
            label = r'$k_T$ = ' + '%.2g' % (k_pre/k) + r' nt/s'
            try:
                localss_population_processing(prefix)
                labels.append(data_ploting(ax_punbound, prefix, label, SD_start+base_position, SD_start+base_position+1, color_rank))
                color_rank += 1
            except:
                continue
        '''
        print('k= inf')
        data = defaultdict(np.float)
        local_structure_collection_data = defaultdict(lambda: defaultdict(np.float))
        prefix = PATH + '/k' + 'inf'
        label = r'$k_T/k_f$ = ' + '0'
        localss_population_processing(prefix)
        color_rank += 1
        labels.append(data_ploting(ax_punbound, prefix, label, SD_start+base_position, SD_start+base_position+1, color_rank))
        '''
        prefix = PATH + '/equilibrium'  # Need to be copied here
        label = 'Equilibrium Control'
        labels.append(data_ploting_equ(ax_punbound, mut, prefix, label, SD_start+base_position, SD_start+base_position+1, color_rank + 1))  # Another format

        # ax_punbound.legend(loc='best')
        adjustText.adjust_text(labels, arrowprops=dict(arrowstyle='-', color='0.7'))
        plt.tick_params(axis="both", which="both", bottom="off", top="off",
                        labelbottom="on", left="off", right="off", labelleft="on")

        fig.savefig(PATH + f'/p_unbound_base[{base_gene_position}]_k_tuning_linear.png',  bbox_inches="tight")
        # plt.show()
exit()
