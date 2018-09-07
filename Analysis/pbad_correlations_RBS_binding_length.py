import math
import numpy as np
import scipy.stats
start_base = 21
end_base = 28
km_start = 4
km_end = 19
k_pre = 1e11
km_interval = 1

import matplotlib
import matplotlib.pyplot as plt
# sphinx_gallery_thumbnail_number = 2

def get_pearsonr(x, y, pval=False):
    R, p = scipy.stats.pearsonr(x, y)
    if np.isfinite(R):
        if pval == True:
            return R, p
        else:
            return R
    else:
        if pval == True:
            return 0, 1
        else:
            return 0


def bootstrap_correlation(x, y, fcn, n=1000):
    z = list(zip(x, y))
    p = 0
    R = np.zeros(n)
    if not x:
        return 0, 0, 1
    for i in range(n):
        s = [z[np.random.randint(len(z))] for j in range(len(z))]
        # s = z
        R[i], p_i = fcn([ss[0] for ss in s], [ss[1] for ss in s], pval=True)
        p += p_i / n
    return np.mean(R), np.std(R), p


if __name__ == '__main__':

    pBAD_upstream = 'ATACCCGTTTTTTGGGCTAACAGGAGGAATTACAT'
    genes = ['adk', 'folA']

    # Load data
    # Define N-term/C-term groups
    Ncterm = 3 * 20
    sequences = {}
    trans_rate_es = 11 - (np.arange(km_start, km_end, km_interval))
    rbs_bind_time_es = list(np.arange(0, 200, 10))
    Nterm_data = np.zeros((len(trans_rate_es), len(rbs_bind_time_es)))

    # Calculate ddG
    for e_k in range(km_start, km_end, km_interval):
        k1 = 1 * 10 ** e_k
        k_T = (k_pre/k1)
        for e_length in np.arange(0, 200, 10):
            # rbs_time = 1 * 10 ** (e_time)
            print(k_T, e_length)
            sequences = {}
            for gene in genes:
                with open('20180110_all_' + gene + '_sequences.txt', 'r') as f:
                    for line in f:
                        if len(line) > 1 and line[0] == '>':
                            seqname = line[1:].strip()
                        elif len(line) > 1:
                            sequences[(gene, seqname)] = line.strip()
            with open('Analysis/table_S1.tsv', 'r') as f:
                data = {}
                for line in f:
                    if len(line) > 1 and line[0] == '#':
                        fields = line.split()[1:]
                    elif len(line) > 1 and len(line.split()) > 1:
                        gene, seqname = line.split()[0], line.split()[1]
                        if (gene, seqname) not in sequences:
                            print("Could not find", gene, seqname)
                            raise Exception
                        data[(gene, seqname)] \
                            = {fields[i]: float(line.split()[i]) if line.split()[i] != 'ND' else np.nan \
                               for i in range(2, len(fields))}

            # Define N-term/C-term groups
            Ncterm = 3 * 20
            group_names = ['mut_nterm', 'mut_ncterm', 'mut_cterm']
            groups = {name: [] for name in group_names}
            for k in sorted(data):
                if k[1] != 'WT':
                    if sequences[k][:Ncterm] != sequences[(k[0], 'WT')][
                                                :Ncterm]:  # and sequences[k][Ncterm:] == sequences[(k[0], 'WT')][Ncterm:]:
                        groups['mut_nterm'].append(k)
                    elif sequences[k][:Ncterm] == sequences[(k[0], 'WT')][:Ncterm] and \
                            sequences[k][Ncterm:] != sequences[(k[0], 'WT')][Ncterm:]:
                        groups['mut_cterm'].append(k)
                    else:
                        groups['mut_ncterm'].append(k)

            for d in data:
                # print(d[0] + '_' + d[1] + '/CG5' + '/k' + '%.2g' % k1 + f'/p_unbound/base{start_base}_{end_base}.dat')
                try:
                    with open(d[0] + '_' + d[1] + '/CG5_pool50' + '/k' + '%.2g' % k1 + f'/p_unbound/base{start_base}_{end_base}.dat', 'r+') as p_unbound_in:
                        mutant_pub_data = {}
                        for line in p_unbound_in:
                            length, pub = map(float, line.split())
                            if int(e_length) == length:
                                # mutant_pub_data[length] = pub
                                data[d]['punbound'] = pub
                except:
                    # print(f'no {d}')
                    continue

            # Calculate correlations
            x_all, y_all = [], []
            '''
            for group in group_names:
                x, y = [], []
                for d in sorted(groups[group]):
                    wt_d = (d[0], 'WT')
                    if data[d]['rel_mrna_abundance'] > 0:
                        x.append(data[d]['ddG'] - data[wt_d]['ddG'])
                        y.append(math.log(data[d]['rel_mrna_abundance']))
                x_all += x
                y += y
                R_mean, R_std, p = bootstrap_correlation(x, y, get_pearsonr)
                print("%s %d %6.3f %6.3f %g" % (k, len(y), R_mean, R_std, p))
            R_mean, R_std, p = bootstrap_correlation(x_all, y_all, get_pearsonr)
            print("all %d %6.3f %6.3f %g" % (len(y), R_mean, R_std, p))
            '''
            x, y = [], []
            for d in sorted(groups['mut_nterm']):
                wt_d = (d[0], 'WT')
                if 'punbound' in data[d] and 'punbound' in data[wt_d] and data[d]['rel_mrna_abundance'] > 0:
                    x.append(math.log(data[d]['punbound']) - math.log(data[wt_d]['punbound']))
                    y.append(math.log(data[d]['rel_mrna_abundance']))
            x_all += x
            y_all += y
            R_mean, R_std, p = bootstrap_correlation(x, y, get_pearsonr)
            print("%s %d %6.3f %6.3f %g" % (d, len(y), R_mean, R_std, p))
            Nterm_data[int((e_k-km_start)/km_interval)][int(e_length/10)] = R_mean



    fig = plt.figure(figsize=(20, 9))
    # colors = [plt.cm.jet(lt) for lt in range(0, 8)]
    fig.add_axes()
    ax = fig.add_subplot(111)
    im = ax.imshow(Nterm_data)

    # We want to show all ticks...
    ax.set_yticks(np.arange(len(trans_rate_es)))
    ax.set_xticks(np.arange(len(rbs_bind_time_es)))
    # ... and label them with the respective list entries
    ax.set_yticklabels(trans_rate_es)
    ax.set_xticklabels(rbs_bind_time_es)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(trans_rate_es)):
        for j in range(len(rbs_bind_time_es)):
            text = ax.text(j, i, '%.2f' % Nterm_data[i, j],
                           ha="center", va="center", color="w")

    ax.set_title("mRNA - p_unbound(effective RBS binding time) correlation")
    fig.tight_layout()
    fig.savefig('Analysis/mRNA-p_unbound(effective RBS binding length)correlation_SD-refined_090718.png')
