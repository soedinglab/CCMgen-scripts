#!/usr/bin/env python

# ===============================================================================
###     This script plots a scatter plot
###     of the apc and ec correction terms for all pairs of residues
###     given a binary raw file defining the MRF model
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

def plot_scatter(apc_mat, ec_mat, plot_file):

    indices_i, indices_j = np.triu_indices(apc_mat.shape[0], k=1)
    apc = apc_mat[indices_i, indices_j]
    ec = ec_mat[indices_i, indices_j]

    text = ["i: {0}<br>j: {1}<br>apc:{2}<br>ec:{3}".format(
        i,j,apc_mat[i, j], ec_mat[i, j])
        for i,j in zip(indices_i, indices_j)]


    scatter_data = go.Scatter(
            x = apc,
            y = ec,
            mode = 'markers',
            marker = dict(color="black"),
            text = text,
            showlegend = False
        )

    diagonal = go.Scatter(
        x=[0, np.max(list(apc) + list(ec))],
        y=[0,np.max(list(apc) + list(ec))],
        mode="lines",
        line=dict(color="darkgrey", width=4, dash="dot"),
        showlegend=False
    )

    pearson_r = pearsonr(apc, ec)

    data=[]
    data.append(diagonal)
    data.append(scatter_data)

    fig = {
        "data": data,
        "layout" : go.Layout(
            hovermode="text",
            font = dict(size=24),
            yaxis = dict(
                title="Entropy Correction",
                exponentformat="e",
                showexponent='All',
                scaleratio=1,
                scaleanchor='x'
            ),
            xaxis = dict(
                title="Average Product Correction",
                exponentformat="e",
                showexponent='All',
                scaleratio=1,
                scaleanchor='y'
            ),
            annotations=go.Annotations([
                go.Annotation(
                    x=0.05,
                    y=0.95,
                    showarrow=False,
                    text='Pearson r = {0}'.format(np.round(pearson_r[0], decimals=3)),
                    font=dict(color="black", size=24),
                    xref='paper',
                    yref='paper'
                )
            ]),
            margin=dict(t=10),
            width="550",
            height="500"
        )
    }


    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)

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
    plot_dir = data_dir +  "/plots/apc_vs_ec/"

    for alignment_file in glob.glob(alignment_dir + "/*aln"):

        protein  = os.path.basename(alignment_file).split(".")[0]

        #PLL
        plot_file  = plot_dir + protein + ".apc_vs_ec.pll.html"
        binary_raw_file = pll_dir + protein + ".braw.gz"
        if os.path.exists(binary_raw_file):
            apc, entropy = compute_correction_terms(alignment_file, binary_raw_file)
            plot_scatter(apc, entropy, plot_file)

        #PCD
        plot_file  = plot_dir + protein + ".apc_vs_ec.pcd.html"
        binary_raw_file = pcd_dir + protein + ".braw.gz"
        if os.path.exists(binary_raw_file):
            apc, entropy = compute_correction_terms(alignment_file, binary_raw_file)
            plot_scatter(apc, entropy, plot_file)


if __name__ == '__main__':
    main()