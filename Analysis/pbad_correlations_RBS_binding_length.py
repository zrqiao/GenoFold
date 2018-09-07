import math
import numpy as np
import scipy.stats
import seaborn as sns; sns.set()
start_base = 21
end_base = 28
km_start = 4
km_end = 19
k_pre = 1e11
km_interval = 1

import matplotlib
import matplotlib.pyplot as plt
# sphinx_gallery_thumbnail_number = 2

def func(x, pos):
    return "{:.2f}".format(x).replace("0.", ".").replace("1.00", "")

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    #for edge, spine in ax.spines.items():
    #    spine.set_visible(False)

    # ax.set_xticks(np.arange(data.shape[1]))
    # ax.set_yticks(np.arange(data.shape[0]))
    # ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    # ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

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
    if len(x) <= 5:
        return float('NaN'), float('NaN'), 1
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
    trans_rate_es = []
    for e_kt in list(11 - (np.arange(km_start, km_end, km_interval))):
        trans_rate_es.append(r'log(Transcription rate)=%g' % e_kt)
    trans_rate_es.append('Equilibrium')
    rbs_bind_time_es = []
    l_rbss = list(np.arange(0, 150, 5))
    for l_rbs in l_rbss: # NOTE: to modify the axes, change here
        rbs_bind_time_es.append(r'Transcript Length=%g' % l_rbs)
    Nterm_data = np.zeros((len(trans_rate_es), len(rbs_bind_time_es)))

    # Calculate ddG
    for e_k in range(km_start, km_end, km_interval):
        k1 = 1 * 10 ** e_k
        k_T = (k_pre/k1)
        for e_length in l_rbss:
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
            Nterm_data[int((e_k-km_start)/km_interval)][int(e_length/5)] = R_mean  # NOTE: Change here

    for e_length in np.arange(0, 150, 5):
        # rbs_time = 1 * 10 ** (e_time)
        print('Equilibrium', e_length)
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
                with open(d[0] + '_' + d[
                    1] + '/CG5_pool50' + '/equilibrium' + f'/p_unbound/base{start_base}_{end_base}.dat',
                          'r+') as p_unbound_in:
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
        Nterm_data[-1][int(e_length / 5)] = R_mean  # NOTE: Change here

    plt.style.use('seaborn')
    fig = plt.figure(figsize=(24, 9))
    # colors = [plt.cm.jet(lt) for lt in range(0, 8)]
    fig.add_axes()
    ax = fig.add_subplot(111)
    # ax = ax.imshow(Nterm_data)
    im, _ = heatmap(Nterm_data, trans_rate_es, rbs_bind_time_es, ax=ax, vmin=-1, vmax=1,
                    cmap="coolwarm", cbarlabel="correlation coeff.")
    #texts = annotate_heatmap(im, valfmt=matplotlib.ticker.FuncFormatter(func))

    # We want to show all ticks...
    #ax.set_yticks(np.arange(len(trans_rate_es)))
    #ax.set_xticks(np.arange(len(rbs_bind_time_es)))
    # ... and label them with the respective list entries
    #ax.set_yticklabels(trans_rate_es)
    #ax.set_xticklabels(rbs_bind_time_es)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=45, ha="right",rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    '''
    for i in range(len(trans_rate_es)):
        for j in range(len(rbs_bind_time_es)):
            text = ax.text(j, i, '%.2f' % Nterm_data[i, j],
                           ha="center", va="center", color="w")
    '''
    ax.set_title("mRNA - "+r"p_unbound"+"(effective RBS binding length) correlation")
    fig.tight_layout()
    fig.savefig('Analysis/mRNA-p_unbound(effective RBS binding length)correlation_SD_090718.png')
