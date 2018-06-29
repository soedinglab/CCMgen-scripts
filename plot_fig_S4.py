#!/usr/bin/env python

# ===============================================================================
###     This script reproduces supplemental Figure 4:
###     It generates a boxplot visualizing the pearson correlation coefficients
###     between the apc and ec correction terms for all pairs of residues
###     over all proteins in the PSICOV dataset.
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import glob
import numpy as np
import ccmpred.raw as raw
from ccmpred import CCMpred
from ccmpred.io import contactmatrix
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from scipy.stats import pearsonr
import sys

def compute_correction_terms(alignment_file, binary_raw_file):

    #initialise ccmpred object
    ccm = CCMpred()

    # specify possible file paths
    ccm.set_alignment_file(alignment_file)
    ccm.set_initraw_file(binary_raw_file)

    # read alignment and remove gapped sequences and positions
    ccm.read_alignment("psicov", 50, 75)

    # compute sequence weights (in order to reduce sampling bias)
    ccm.compute_sequence_weights("simple", 0.8)

    # compute amino acid counts and frequencies adding pseudo counts for non-observed amino acids
    ccm.compute_frequencies("uniform_pseudocounts", 1,  1)

    #read in binary raw file
    ccm.intialise_potentials()

    #compute apc
    ccm.recenter_potentials()
    cmat = contactmatrix.frobenius_score(ccm.x_pair)
    mean = np.mean(cmat, axis=0)
    apc_mat = mean[:, np.newaxis] * mean[np.newaxis, :] / np.mean(cmat)

    #compute entropy correction
    single_freq = ccm.pseudocounts.freqs[0]
    nr_states = 20
    log = np.log2
    scaling_factor, mat_corrected = contactmatrix.compute_local_correction(
        single_freq, ccm.x_pair, ccm.neff, 1,
        squared=False, entropy=True, nr_states=nr_states, log=log
    )
    entropy_correction_mat = cmat - mat_corrected

    return apc_mat, entropy_correction_mat

def plot_boxplot_correlation(pearson_r_pll, pearson_r_pcd, plot_file):


    box_pearson_pll = go.Box(
        y=pearson_r_pll,
        name = "pseudo-likelihood",
        showlegend=False,
        boxmean=False,
        boxpoints='Outliers'
    )

    box_pearson_pcd = go.Box(
        y=pearson_r_pcd,
        name="persistent contrastive divergence",
        showlegend=False,
        boxmean=False,
        boxpoints='Outliers'
    )


    plot = {
    "data": [box_pearson_pll, box_pearson_pcd],
    "layout" : go.Layout(
        title = "Correlation between APC and Entropy Correction",
        font = dict(size=24),
        margin=dict(t=50),
        yaxis=dict(range=[0,1], title="Pearson correlation"),
        width="900",
        height="450"
        )
    }

    plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)

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

    pll_dir = data_dir + "/predictions_pll/"
    pcd_dir = data_dir + "/predictions_pcd/"
    alignment_dir = data_dir  + "/aln/"
    plot_dir = data_dir + "/plots/supplement/"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    if not os.path.exists(pll_dir) or not os.path.exists(pcd_dir):
        print("You first need to learn Markov Random Field models in {0} and {1}".format(
            pll_dir, pcd_dir))
        sys.exit(1)


    pearson_r_list_pll = []
    pearson_r_list_pcd = []
    for alignment_file in glob.glob(alignment_dir + "/*aln"):

        protein  = os.path.basename(alignment_file).split(".")[0]

        #PLL
        binary_raw_file = pll_dir + protein + ".braw.gz"
        if os.path.exists(binary_raw_file):

            try:
                apc, entropy = compute_correction_terms(alignment_file, binary_raw_file)
                indices_i, indices_j = np.triu_indices(apc.shape[0], k=1)
                # compute pearson correlation coefficient
                pearson_r_list_pll.append(
                    pearsonr(apc[indices_i, indices_j], entropy[indices_i, indices_j])[0])
            except:
                print("Unexpected error:", sys.exc_info()[0])



        #PCD
        binary_raw_file = pcd_dir + protein + ".braw.gz"
        if os.path.exists(binary_raw_file):
            try:
                apc, entropy = compute_correction_terms(alignment_file, binary_raw_file)
                indices_i, indices_j = np.triu_indices(apc.shape[0], k=1)

                # compute pearson correlation coefficient
                pearson_r_list_pcd.append(
                    pearsonr(apc[indices_i, indices_j], entropy[indices_i, indices_j])[0])
            except:
                print("Unexpected error:", sys.exc_info()[0])



    plot_file = plot_dir + "/fig_S4.html"
    plot_boxplot_correlation(pearson_r_list_pll, pearson_r_list_pcd, plot_file)

if __name__ == '__main__':
    main()