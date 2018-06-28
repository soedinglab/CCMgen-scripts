#!/usr/bin/env python

# ===============================================================================
###     This script plots Figure 1c of the CCMgen paper
###     Contact Prediction Benchmark (mean precision vs best predictions)
###     for contact predictions from MRF models inferred with
###     - pseudo-likelihood maximiation
###     - persistent contrastive divergence
# ===============================================================================

### load libraries
import argparse
import os
from benchmark import Benchmark
import copy
from plotly.offline import plot as plotly_plot

def plot_pll_vs_pcd_benchmark_figure(plot, plot_dir, height=400, width=600):

    data = []

    #first add a legend for the methods
    for trace in plot['data']:
        #adding a legend for the methods
        if "persistent contrastive diverence APC" in trace['name']:
            trace_for_legend = copy.copy(trace)
            data.append(trace_for_legend)
            data[-1]['name'] = 'persistent contrastive divergence'
            data[-1]['legendgroup'] = 'method'
            data[-1]['line']['color'] = 'black'
            data[-1]['showlegend'] = True

        if "pseudo-likelihood APC" in trace['name']:
            trace_for_legend = copy.copy(trace)
            data.append(trace_for_legend)
            data[-1]['name'] = 'pseudo-likelihood'
            data[-1]['legendgroup'] = 'method'
            data[-1]['line']['color'] = 'black'
            data[-1]['line']['dash'] = 'dot'
            data[-1]['showlegend'] = True

    #add all other traces and legend for correction
    for trace in plot['data']:
        trace['legendgroup']='correction'

        #pseudo-likelihood will be dotted
        if "pseudo-likelihood" in trace['name']:
            trace['line']['dash'] = 'dot'
            trace['showlegend'] = False

        #formatting legend for corrections
        if "APC" in trace['name']:
            trace['line']['color'] = 'blue'
            trace['name'] = "APC"

        if "raw" in trace['name']:
            trace['line']['color'] = 'red'
            trace['name'] = "no APC"

        data.append(trace)

    #replace traces in plot
    plot['data'] = data

    plot['layout']['legend'] = {
        "orientation": "v",
        "x":1.01,
        "y":1.0
    }
    plot['layout']['height'] = height
    plot['layout']['width'] = width

    plot_file = plot_dir+"/"+"fig_1c.html"
    print("Write benchmark plot to {0}.".format(plot_file))
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
    plot_dir = data_dir+"/plots/benchmarks/"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    # create benchmark object
    b = Benchmark(data_dir+"/pdb/")


    #specify methods to benchmark
    b.add_method("pseudo-likelihood APC", data_dir +"/predictions_pll/", "apc.mat")
    b.add_method("pseudo-likelihood raw", data_dir +"/predictions_pll/", "raw.mat")
    b.add_method("persistent contrastive diverence APC", data_dir +"/predictions_pcd/", "apc.mat")
    b.add_method("persistent contrastive diverence raw", data_dir +"/predictions_pcd/", "raw.mat")

    #add constraint that all MRF optimizations have exist status 0
    b.add_constraint("opt_code", 0, "greater_equal")

    #compute the precision of predictions
    b.compute_evaluation_statistics(seqsep=6, contact_thr=8, noncontact_thr=8)

    #generate a benchmark plot
    plot = b.plot_precision_vs_rank()

    #format that benchmark plot to resemble the one in Fig 1C
    plot_pll_vs_pcd_benchmark_figure(plot, plot_dir, height=500, width=1000)




if __name__ == '__main__':
    main()