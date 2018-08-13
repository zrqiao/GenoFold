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

    #   NOTE: Initiation [create active population]

    # NOTE: Initiate transcription

    # print('Population size: '+str(active_species_pool.size))

    # Start IO
    fig = plt.figure(figsize=(8, 6))
    fig.add_axes()
    # gs = gridspec.GridSpec(1, 1, height_ratios=[1])
    NUM_COLORS = 16
    #ax_energy = fig.add_subplot(gs[0, 0])
    #ax_energy.set_title('Free Energy')
    # ax_energy.set_xlabel('Subsequence Length', fontsize=12.5)
    #ax_energy.set_ylabel('Free Energy', fontsize=12.5)
    cm = plt.get_cmap('gist_rainbow')
    ax_pbound = fig.add_subplot(111)
    ax_pbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

    # ax_pbound.set_title(f'Average p_unbound for base[-9](G) {clargs.sequence}')
    ax_pbound.set_title(f'Average p_unbound of SD sequence {clargs.sequence}')
    ax_pbound.set_xlabel('Transcription time', fontsize=12.5)
    ax_pbound.set_ylabel(r'$p_{unbound}$', fontsize=12.5)
    # ax_pbound.set_yscale('log')
    # ax_pbound.set_xscale('log')
    # ax_pbound.set_ylim(1e-5, 1.5)
    # ax_pbound.set_ylim(0.0, 1.1)

    for e_k in range(1, 16, 1):
        k = 1*10**e_k
        print(f'k= {k}')
        data = defaultdict(np.float)
        f = open(clargs.sequence + '_p_unbound_%e' % k + '.dat', 'w')
        with open(clargs.sequence + '_k' + '%e' % k + '.dat', 'r+') as folding_input:

            sss = [(x.split()[0], np.float(x.split()[1]))
                   for x in folding_input.readlines()]

            for ss in sss:
                if ss[0].startswith('#'):
                    time = ss[1]
                else:
                    # print(ss)
                    if len(ss[0]) >= SD_end:
                        # time = (len(ss[0])-L_init) * transcription_time
                        # print('N=' + str(N) + ' nt')
                        # f.write('#N=' + str(N) + ' nt')
                        # sim = np.zeros(len(structure))
                        SD_ss = ss[0][SD_start:SD_end]

                        data[time] += ss[1]  # * similar(SD_ss, '.......')
                        # print(data[time])
                        # for i in range(len(first_sequences)):
                        #     temp_cmp = []
                        #     for cmpseq in first_sequences:
                        #         temp_cmp.append(similar(cmpseq, first_sequences[i]))
                        #     sim[i + N] = np.average(np.array(temp_cmp))
                        # print(sim)
            data_p = np.array([list(data.keys()), list(data.values())])
            # data_p.transpose()
            # folding_input.close()

        ax_pbound.plot(data_p[0], data_p[1], label='k=' + '%.2g' % k + r'$s^{-1}$')
        for d in data.items():
            f.write(f'{d[0]}  {d[1]}\n')
        f.close()
    ax_pbound.plot(data_p[0], np.repeat(np.average(equi_p_unbound), len(data_p[0])), label='equilibrium')
    ax_pbound.legend(loc='best')
    # fig.tight_layout()
    plt.show()

    fig.savefig(clargs.sequence + f'_p_unbound_SD_k_tuning.eps')

    for base_position in range(0, 7):
        base_gene_position = base_position-14
        fig = plt.figure(figsize=(8, 6))
        fig.add_axes()
        # gs = gridspec.GridSpec(1, 1, height_ratios=[1])
        # NUM_COLORS = 13
        # ax_energy = fig.add_subplot(gs[0, 0])
        # ax_energy.set_title('Free Energy')
        # ax_energy.set_xlabel('Subsequence Length', fontsize=12.5)
        # ax_energy.set_ylabel('Free Energy', fontsize=12.5)
        ax_pbound = fig.add_subplot(111)
        ax_pbound.set_color_cycle([cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])

        # ax_pbound.set_title(f'Average p_unbound for base[-9](G) {clargs.sequence}')
        ax_pbound.set_title(f'p_unbound of base [{base_gene_position}] {clargs.sequence}')
        ax_pbound.set_xlabel('Transcription time', fontsize=12.5)
        ax_pbound.set_ylabel(r'$p_{unbound}$', fontsize=12.5)
        ax_pbound.set_yscale('log')
        # ax_pbound.set_xscale('log')
        ax_pbound.set_ylim(1e-4, 1.5)

        for e_k in range(1, 16, 1):
            k = 1 * 10 ** e_k
            print(f'k= {k}')
            data = defaultdict(np.float)
            f = open(clargs.sequence + '_p_unbound_%e' % k + '.dat', 'w')
            with open(clargs.sequence + '_k' + '%e' % k + '.dat', 'r+') as folding_input:
                sss = [(x.split()[0].rstrip('\n'), np.float(x.split()[1]))
                       for x in folding_input.readlines()]

                for ss in sss:
                    if ss[0].startswith('#'):
                        time = ss[1]
                    else:
                        # print(time)
                        if len(ss[0]) >= SD_end:
                            # time = (len(ss[0])-L_init) * transcription_time
                            # print('N=' + str(N) + ' nt')
                            # f.write('#N=' + str(N) + ' nt')
                            # sim = np.zeros(len(structure))
                            SD_ss = ss[0][SD_start:SD_end]
                            data[time] += ss[1] * similar(SD_ss[base_position], '.')
                            # print(data[time])
                            # for i in range(len(first_sequences)):
                            #     temp_cmp = []
                            #     for cmpseq in first_sequences:
                            #         temp_cmp.append(similar(cmpseq, first_sequences[i]))
                            #     sim[i + N] = np.average(np.array(temp_cmp))
                            # print(sim)
                data_p = np.array([list(data.keys()), list(data.values())])
                # data_p.transpose()
            # print(data_p)
            ax_pbound.plot(data_p[0], data_p[1], label='k=' + '%.2g' % k + r'$s^{-1}$')
            for d in data.items():
                f.write(f'{d[0]}  {d[1]}\n')
            f.close()
        ax_pbound.plot(data_p[0], np.repeat(equi_p_unbound[base_position],
                                            len(data_p[0])), label='equilibrium')
        ax_pbound.legend(loc='best')
        # fig.tight_layout()
        plt.show()

        fig.savefig(clargs.sequence + f'_p_unbound_base[{base_gene_position}]_k_tuning.eps')
exit()
