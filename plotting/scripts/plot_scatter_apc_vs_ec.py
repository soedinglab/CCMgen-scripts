#!/usr/bin/env python

# ===============================================================================
###     This script plots a scatter plot
###     of the apc and ec correction terms for all pairs of residues
###     given a binary raw file defining the MRF model
# ===============================================================================

### load libraries ===============================================================================
import argparse
import os
import numpy as np
import ccmpred.raw as raw
import plotting.utils.io_utils as io
from ccmpred import CCMpred
from ccmpred.io import contactmatrix
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from scipy.stats import pearsonr


def plot_scatter(apc, ec, text, plot_file):

    scatter_data = go.Scatter(
            x= apc,
            y= ec,
            mode = 'markers',
            marker=dict(color="black"),
            text = text,
            showlegend=False
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
    Parse command line arguments

    :return:
    """

    parser = argparse.ArgumentParser(description='Plot scatter plot of apc vs ec correction terms.')
    parser.add_argument("binary_raw_file",   type=str, help="path to binary raw file")
    parser.add_argument("alignment_file",   type=str, help="path to alignment file")
    parser.add_argument("plot_file",         type=str, help="path to output plot")

    args = parser.parse_args()

    return args

def main():


    args = parse_args()

    binary_raw_file = args.binary_raw_file
    alignment_file = args.alignment_file
    plot_file  = args.plot_file

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
    apc_term = mean[:, np.newaxis] * mean[np.newaxis, :] / np.mean(cmat)

    #compute entropy correction
    single_freq = ccm.pseudocounts.freqs[0]
    nr_states = 20
    log = np.log2
    scaling_factor, mat_corrected = contactmatrix.compute_local_correction(
        single_freq, ccm.x_pair, ccm.neff, 1,
        squared=False, entropy=True, nr_states=nr_states, log=log
    )
    entropy_correction = cmat - mat_corrected

    #generate the scatter plot
    indices_i, indices_j = np.triu_indices(ccm.L, k=1)
    plot_scatter(
        apc_term[indices_i, indices_j],
        entropy_correction[indices_i, indices_j],
        ["i: " + str(i) + "<br>j: " + str(j) for i,j in zip(indices_i, indices_j)],
        plot_file)






if __name__ == '__main__':
    main()