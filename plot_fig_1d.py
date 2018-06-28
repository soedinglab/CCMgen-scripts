#!/usr/bin/env python

# ===============================================================================
###     This script plots a boxplot of the runtimes in minutes
###     for CCMpredPy runs from different methods/settings
# ===============================================================================

### load libraries
import argparse
import sys
import os
import glob
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import ccmpred.io.contactmatrix

def parse_args():
    """
    parse command line arguments
    :return:
    """

    parser = argparse.ArgumentParser(description='Plot CCMgen paper Figure 1C.')
    parser.add_argument("data_dir", type=str, help="path to psicov data working directory")

    args = parser.parse_args()

    return args

def plot_runtime(plot_data, plot_file):

    data = []
    for method,runtimes in plot_data.items():

        box = go.Box(
            y=runtimes,
            boxmean=True,
            boxpoints='Outliers',
            name=method,
            marker=dict(opacity=1),
            hoverinfo='all',
            orientation='v',
            showlegend=False
        )

        data.append(box)


    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                title="runtime in min",
                type='log',
                exponentformat='none',
                showexponent='none',
                tickmode="array",
                tickvals=[1, 10, 100, 500, 1000, 5000, 10000],
                ticktext=[1, 10, 100, 500, 1000, 5000, 10000]
            ),
            font=dict(size=18)
        )
    }

    plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)

def main():

    #parse command line arguments
    args = parse_args()

    data_dir = args.data_dir
    plot_dir = data_dir + "/plots/"
    mat_dirs = [data_dir + "/predictions_pll/", data_dir + "/predictions_pcd/"]
    methods = ["pseudo-likelihood", "persistent contrastive divergence"]


    if not os.path.exists(plot_dir):
        print("Plot dir {0} does not exitst!".format(plot_dir))
        sys.exit(1)


    #prepare dictionary for plotting
    plot_data = {}

    #iterate over all contact matrix files for all methods
    for id, mat_dir in enumerate(mat_dirs):
        mat_files = glob.glob(mat_dir + "/*.apc.mat")
        plot_data[methods[id]] = []

        for mat_file in mat_files:
            mat, meta = ccmpred.io.contactmatrix.read_matrix(mat_file)
            runtime =ccmpred.io.contactmatrix.find_dict_key("runtime", meta)
            plot_data[methods[id]].append(runtime)


    #plot the distribution of runtimes for each method as boxplot
    plot_file = plot_dir + "fig_1d.html"
    plot_runtime(plot_data, plot_file)



if __name__ == '__main__':
    main()
