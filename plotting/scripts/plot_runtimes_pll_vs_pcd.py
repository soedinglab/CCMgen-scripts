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
import plotting.utils.io_utils as io

def parse_args():
    """
    Parse command line arguments
    :return:
    """

    parser = argparse.ArgumentParser(description='plot benchmark for specified eval files and scores.')
    parser.add_argument("plot_dir",  type=str, help="path to print plot files")
    parser.add_argument("mat_dirs",  type=str, help="comma separated string of paths to folders with *mat files")
    parser.add_argument("method",    type=str, help="comma separated string of method names")

    args = parser.parse_args()

    return args

def plot_runtime(plot_data, methods, plot_dir):

    plot_name = plot_dir+"/runtime_boxplot"
    for method in methods:
        plot_name+="_"+str(method)
    plot_name+=".html"


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

    plotly_plot(plot, filename=plot_name, auto_open=False, show_link=False)



def main():

    #parse command line arguments
    args = parse_args()

    plot_dir        = args.plot_dir
    mat_dirs        = args.mat_dirs.split(",")
    methods         = args.method.split(",")

    if len(mat_dirs) != len(methods):
        print("Number of specified method names does not match number of specified mat folders!")
        sys.exit(1)


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
            meta = io.read_json_from_mat(mat_file)
            runtime =io.find_dict_key("runtime", meta)
            plot_data[methods[id]].append(runtime)



    #plot the distribution of runtimes for each method as boxplot
    plot_runtime(plot_data, methods, plot_dir)



if __name__ == '__main__':
    main()
