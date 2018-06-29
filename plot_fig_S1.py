#!/usr/bin/env python

# ===============================================================================
###     This script reproduces supplemental Figure 1:
###     It generates boxplots visualizing the pearson correlation coefficients
###     between the alignment statistics from the original Pfam alignment and
###     MCMC samples drawn from either a pseudo-likelihood MRF model or a MRF learned with PCD
###     over all proteins in the PSICOV dataset.
# ===============================================================================

### load libraries ===============================================================================
import argparse
import sys
import os
import glob
import numpy as np
from scipy.stats import pearsonr

import ccmpred.io.alignment
import ccmpred.gaps
import ccmpred.pseudocounts
import ccmpred.weighting

import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot

def plot_boxplot_correlation_alignment_statistics_pll_vs_pcd(data_dict, plot_file):

    data = []

    data.append(
        go.Box(
            y=data_dict['pseudo-likelihood']['x'],
            x=data_dict['pseudo-likelihood']['y'],
            boxpoints='outliers',
            name="pseudo-likelihood",
            hoverinfo='all',
            orientation="v",
            showlegend=True
        )
    )

    data.append(
        go.Box(
            y=data_dict['contrastive divergence']['x'],
            x=data_dict['contrastive divergence']['y'],
            boxpoints='outliers',
            name="persistent contrastive divergence",
            hoverinfo='all',
            orientation="v",
            showlegend=True
        )
    )

    layout=go.Layout(
        #title="Pearson Correlation Coefficients<br>between Original and Sampled Alignment Statistics",
        title="",
        margin=dict(t=10),
        legend=dict(orientation="h",
                    xanchor="center", x=0.5, y=1.2),
        yaxis=dict(title="Pearson's r", range=[0,1]),
        font=dict(size=18),
        boxmode='group'
    )

    fig = go.Figure(data=data, layout=layout)

    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)

def get_freq(alignment):

    # compute sequence weights for observed sequences
    weights = ccmpred.weighting.weights_simple(alignment, 0.8)

    # compute observed amino acid frequencies
    pseudocounts = ccmpred.pseudocounts.PseudoCounts(alignment, weights)
    pseudocounts.calculate_frequencies(
        'uniform_pseudocounts', 1, 1, remove_gaps=False
    )
    single_freq, pairwise_freq = pseudocounts.freqs

    # degap the frequencies (ignore gap frequencies)
    single_freq = pseudocounts.degap(single_freq, False)
    pairwise_freq = pseudocounts.degap(pairwise_freq, False)

    return single_freq, pairwise_freq

def parse_args():
    """
    parse command line arguments
    :return:
    """

    parser = argparse.ArgumentParser(description='Plot CCMgen paper Figure 1C.')
    parser.add_argument("data_dir", type=str, help="path to psicov data working directory")

    args = parser.parse_args()

    return args

def main():

    #parse command line arguments
    args = parse_args()

    data_dir = args.data_dir

    plot_dir = data_dir + "/plots/supplement/"
    alignment_dir = data_dir + "/aln/"
    samples_pll_dir = data_dir + "/samples_pll/"
    samples_pcd_dir = data_dir + "/samples_pcd/"
    max_gap_pos = 50

    if not os.path.exists(samples_pll_dir) or not os.path.exists(samples_pcd_dir):
        print("You first need to generate MCMC samples from Markov Random Field models in {0} and {1}".format(
            samples_pll_dir, samples_pcd_dir))
        sys.exit(1)

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    data_dict = {
        'pseudo-likelihood': {
            'x': [],
            'y': []
        },
        'contrastive divergence': {
            'x': [],
            'y': []
        }
    }


    for alignment_file in glob.glob(alignment_dir + "/*aln"):

        protein  = os.path.basename(alignment_file).split(".")[0]

        sampled_pll = samples_pll_dir + "/" + protein + ".mcmc.aln"
        sampled_pcd = samples_pcd_dir + "/" + protein + ".mcmc.aln"

        if not os.path.exists(sampled_pll) or not os.path.exists(sampled_pcd):
            continue

        print("compute correlation of alignment statistics for {0}...".format(protein))

        #read alignment files
        alignment = ccmpred.io.alignment.read_msa_psicov(alignment_file)
        alignment_sampled_pll = ccmpred.io.alignment.read_msa_psicov(sampled_pll)
        alignment_sampled_pcd = ccmpred.io.alignment.read_msa_psicov(sampled_pcd)

        #remove gappy positions
        alignment, gapped_positions = ccmpred.gaps.remove_gapped_positions(alignment, max_gap_pos)
        non_gapped_positions = [i for i in range(alignment_sampled_pll.shape[1]) if i not in gapped_positions]
        alignment_sampled_pll = np.ascontiguousarray(alignment_sampled_pll[:, non_gapped_positions])
        alignment_sampled_pcd = np.ascontiguousarray(alignment_sampled_pcd[:, non_gapped_positions])

        #get amino acid frequencies
        freq_single, freq_pair = get_freq(alignment)
        freq_single_pll, freq_pair_pll = get_freq(alignment_sampled_pll)
        freq_single_pcd, freq_pair_pcd = get_freq(alignment_sampled_pcd)


        #reshape and compute covariances
        L = alignment.shape[1]
        indices_i, indices_j = np.triu_indices(L, k=1)

        single = freq_single.flatten().tolist()
        single_pll = freq_single_pll.flatten().tolist()
        single_pcd = freq_single_pcd.flatten().tolist()

        pair = freq_pair[indices_i, indices_j, :, :].flatten().tolist()
        pair_pll = freq_pair_pll[indices_i, indices_j, :, :].flatten().tolist()
        pair_pcd = freq_pair_pcd[indices_i, indices_j, :, :].flatten().tolist()

        cov = [freq_pair[i, j, a, b] - (freq_single[i, a] * freq_single[j, b])
                        for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        cov_pll = [freq_pair_pll[i, j, a, b] - (freq_single_pll[i, a] * freq_single_pll[j, b])
                        for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]
        cov_pcd = [freq_pair_pcd[i, j, a, b] - (freq_single_pcd[i, a] * freq_single_pcd[j, b])
                        for i in range(L - 1) for j in range(i + 1, L) for a in range(20) for b in range(20)]


        #compute pearson correlation
        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(single, single_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('single site amino<br>acid frequencies')

        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(pair, pair_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('pairwise amino<br>acid frequencies')

        data_dict['pseudo-likelihood']['x'].append(np.corrcoef(cov, cov_pll)[0, 1])
        data_dict['pseudo-likelihood']['y'].append('Covariances')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(single, single_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('single site amino<br>acid frequencies')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(pair, pair_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('pairwise amino<br>acid frequencies')

        data_dict['contrastive divergence']['x'].append(np.corrcoef(cov, cov_pcd)[0, 1])
        data_dict['contrastive divergence']['y'].append('Covariances')



    #plot boxplot
    plot_file = plot_dir + "/fig_S1.html"
    plot_boxplot_correlation_alignment_statistics_pll_vs_pcd(data_dict, plot_file)

if __name__ == '__main__':
    main()