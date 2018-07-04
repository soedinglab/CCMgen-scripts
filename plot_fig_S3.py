#!/usr/bin/env python

# ===============================================================================
###     This script plots Supplemental Figures 3a and 3b of the CCMgen paper:
###     a) Scatter plots of coupling scores learned with pseudo-likelihood
###         maximization and persistent contrastive divergence or every
###         protein in the PSICOV dat aset
###     b) box plots visualizing the distribution of correlation statistics
###         between coupling scores computed form both methods over all
###         proteins in the data set
# ===============================================================================

### load libraries
import argparse
import os
from plotly.offline import plot as plotly_plot
import plotly.graph_objs as go
import glob
from ccmpred.io import contactmatrix
import numpy as np
from scipy.stats import ks_2samp, spearmanr, kendalltau, pearsonr, linregress
import pandas as pd

def plot_scatter_comparison(mat_pll, mat_pcd, plot_file, qqplot=False):

    L = mat_pll.shape[0]
    indices_i, indices_j = np.triu_indices(L, k=1)

    #get contact score vectors
    score_pll = mat_pll[indices_i, indices_j]
    score_pcd = mat_pcd[indices_i, indices_j]

    #compute linear regression
    lin_reg_x = list(np.arange(
        np.min([np.min(score_pll),np.min(score_pcd)]),
        np.max([np.max(score_pll),np.max(score_pcd)]),
        0.05))
    slope, intercept, rvalue, pvalue, stderr = linregress(score_pcd, score_pll)
    lin_reg_y = [intercept + slope * x for x in lin_reg_x]

    #hover text in interactive html file
    text = ["i: " + str(i+1) + "<br>j: " + str(j+1) for i,j in zip(indices_i, indices_j)]


    data=[]

    #plot diagonal at bottom
    data.append(
        go.Scattergl(
            x=[np.min([np.min(score_pll), np.min(score_pcd)]), np.min([np.max(score_pll), np.max(score_pcd)])],
            y=[np.min([np.min(score_pll), np.min(score_pcd)]), np.min([np.max(score_pll), np.max(score_pcd)])],
            mode='lines',
            line=dict(color='lightgrey',
                      width=3,
                      dash='dash'
                      ),
            showlegend=False
        )
    )

    # plot scatter in blue
    data.append(
        go.Scattergl(
            x= score_pcd,
            y= score_pll,
            text = text,
            mode = 'markers',
            marker=dict(
                opacity=1,
                color="rgb(31,120,180)"
            ),
            hoverinfo="x+y+text",
            showlegend=False
        )
    )

    # QQ Plot will be printed on top in orange
    if qqplot:
        index_sorted_i = np.argsort(score_pll)
        index_sorted_j = np.argsort(score_pcd)

        text_sorted = ["i: " + str(i + 1) + "<br>j: " + str(j + 1) for i, j in zip(
            indices_i[index_sorted_i],
            indices_j[index_sorted_j]
        )]

        data.append(
            go.Scattergl(
                x=sorted(score_pcd),
                y=sorted(score_pll),
                text=text_sorted,
                mode='markers',
                marker=dict(
                    color="rgb(255,127,0)"),
                hoverinfo="x+y+text",
                showlegend=False
            )
        )

    # plot linear regression fit as black line
    data.append(
        go.Scatter(
            x=lin_reg_x,
            y=lin_reg_y,
            mode='lines',
            opacity=1,
            line=dict(
                color="rgb(0,0,0)",
                width=4),
            showlegend=False
        )
    )



    #add regression fit function
    annotation = dict(
        x=np.percentile(lin_reg_x, 95),
        y=np.percentile(lin_reg_y, 95),
        showarrow=True,
        ax=30,
        ay=50,
        arrowcolor='black',
        arrowside="start",
        text='y = {0} + {1}x'.format(np.round(intercept, decimals=3),np.round(slope, decimals=3)),
        font=dict(color="black", size=18),
        xref='x1',
        yref='y1'
    )


    plot = {
        "data": data,
        "layout" : go.Layout(
            annotations=[annotation],
            title = "",
            margin=dict(t=10),
            font=dict(size=18),
            yaxis1 = dict(
                title="pseudo-likelihood maximization",
                exponentformat="e",
                showexponent='all',
                scaleratio=1.0,
                scaleanchor='x',
            ),
            xaxis1 = dict(
                title="persistent contrastive divergence",
                exponentformat="e",
                showexponent='all',
                scaleratio=1.0,
                scaleanchor='y'
            ),
            width=800,
            height=800
        )
    }

    plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)

def plot_boxplot_correlation(stats_dict, keys_list, plot_file):

    df = pd.DataFrame(stats_dict)
    df = df.transpose()

    df['Pearson r'] = [x for x,y in df['pearson'].tolist()]
    df['Pearson pvalue'] = [y for x,y in df['pearson'].tolist()]
    df['Spearman rho'] = [x for x,y in df['spearmanrho'].tolist()]
    df['Spearman pvalue'] = [y for x,y in df['spearmanrho'].tolist()]
    df['Kendalls tau'] = [x for x,y in df['kendalltau'].tolist()]
    df['Kendalls pvalue'] = [y for x,y in df['kendalltau'].tolist()]
    df['kolmogorov-smirnov pvalue'] = [y for x,y in df['kolmogorov-smirnov'].tolist()]
    df['kolmogorov-smirnov'] = [x for x,y in df['kolmogorov-smirnov'].tolist()]
    df['linear fit slope'] = [slope for slope, intercept, rvalue, pvalue, stderr in df['linreg'].tolist()]
    df['linear fit intercept'] = [intercept for slope, intercept, rvalue, pvalue, stderr in df['linreg'].tolist()]

    df['protein'] = df.index

    data = []
    for key in keys_list:
        data.append(
            go.Box(
                y=df[key],
                name = key,
                text=df['protein'],
                showlegend=False,
                boxmean=False,
                boxpoints='outliers'
            )
        )

    plot = {
        "data": data,
        "layout": go.Layout(
            margin=dict(t=10),
            font=dict(size=18),
            yaxis1=dict(
                title="statistics value",
                exponentformat="e",
                showexponent='all',
                range=[0,1]
            ),
            width=800,
            height=500
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


    args = parse_args()

    data_dir = args.data_dir
    plot_dir = data_dir+"/plots/supplement/"
    pll_dir = data_dir+"/predictions_pll/"
    pcd_dir = data_dir + "/predictions_pcd/"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    if not os.path.exists(plot_dir + "/pll_vs_pcd_apc_score_comparison/"):
        os.makedirs(plot_dir + "/pll_vs_pcd_apc_score_comparison/")


    stats_dict={}
    for mat_file_pll in glob.glob(pll_dir + "/*.apc.mat"):

        protein = os.path.basename(mat_file_pll).split(".")[0]
        mat_file_pcd = pcd_dir + "/" + protein + ".apc.mat"

        if not os.path.exists(mat_file_pcd):
            continue

        print("Computing statistics for protein {0}...".format(protein))

        #read contact matrices (corrected with APC)
        mat_pll, meta_pll = contactmatrix.read_matrix(mat_file_pll)
        mat_pcd, meta_pcd = contactmatrix.read_matrix(mat_file_pcd)

        plot_file = plot_dir +  '/pll_vs_pcd_apc_score_comparison/'  + protein  + "_scatter_pll_vs_pcd.html"
        plot_scatter_comparison(mat_pll, mat_pcd, plot_file, qqplot=True)

        #compute correlation statistics
        L = mat_pll.shape[0]
        scores_pll = mat_pll[np.triu_indices(L, k=1)]
        scores_pcd = mat_pcd[np.triu_indices(L, k=1)]
        stats_dict[protein] = {
            "pearson": pearsonr(scores_pcd, scores_pll),
            "kolmogorov-smirnov":  ks_2samp(scores_pcd, scores_pll),
            "spearmanrho": spearmanr(scores_pcd, scores_pll),
            "kendalltau": kendalltau(scores_pcd, scores_pll),
            "linreg": linregress(scores_pcd, scores_pll)
        }

    plot_file = plot_dir +  '/'  + "fig_S3b.html"
    plot_boxplot_correlation(stats_dict, ["Pearson r", "Spearman rho", "Kendalls tau", "linear fit slope"], plot_file)





if __name__ == '__main__':
    main()