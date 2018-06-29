#!/usr/bin/env python

# ===============================================================================
###     This script plots Figure 6 of the CCMgen paper
###     Contact Prediction Benchmark (mean precision vs best predictions)
###     for contact predictions from MRF models inferred with
###     - persistent contrastive divergence
###     from synthetic alignments generated with CCMgen
# ===============================================================================

### load libraries
import argparse
import os
import numpy as np
from benchmark import Benchmark
import copy
from plotly.offline import plot as plotly_plot
import plotly.graph_objs as go

def plot_ccmgen_benchmark_figure(fig, title, plot_file, height=350, width=400):

    fig['layout']['font']['size'] =18
    fig['layout']['hovermode']='closest'
    fig['layout']['title']=title
    fig['layout']['margin']['b']=45
    fig['layout']['margin']['t']=50
    fig['layout']['legend']={
        'orientation':"v",
        'x':0.65, 'y': 1.0
    }
    fig['layout']['yaxis']={
        'title': "mean precision over proteins",
        'range' : [0,1]
    }
    fig['layout']['height'] = height
    fig['layout']['width'] = width

    plotly_plot(fig, filename=plot_file, auto_open=False, show_link=False)

def plot_ccmgen_noise_quant_figure(benchmark_plot_star, benchmark_plot_binary, plot_file, height=350, width=400):


    ### define entropy noise for star topology
    precision_noapc_star = []
    precision_ec_star = []
    x = []
    for trace in benchmark_plot_star['data']:
        if 'no APC' in trace['name']:
            precision_noapc_star = trace['y']
        if 'EC' in trace['name']:
            precision_ec_star = trace['y']
        x = trace['x']

    entropy_noise_star = np.array(precision_ec_star) - np.array(precision_noapc_star)
    entropy_noise_star_trace = go.Scatter(
        x = x,
        y = entropy_noise_star,
        name="entropy noise star",
        line=dict(width=4)
    )

    ### define entropy noise for binary topology
    precision_noapc_binary = []
    precision_ec_binary = []
    for trace in benchmark_plot_binary['data']:
        if 'no APC' in trace['name']:
            precision_noapc_binary = trace['y']
        if 'EC' in trace['name']:
            precision_ec_binary = trace['y']

    entropy_noise_binary = np.array(precision_ec_binary) - np.array(precision_noapc_binary)
    entropy_noise_binary_trace = go.Scatter(
        x = x,
        y = entropy_noise_binary,
        name="entropy noise binary",
        line=dict(width=4)
    )

    ### define phylogenetic noise

    phylogenetic_noise = np.array(precision_ec_star) - np.array(precision_ec_binary)
    phylogenetic_noise_trace = go.Scatter(
        x = x,
        y = phylogenetic_noise,
        name="phylogenetic noise",
        line=dict(width=4)
    )

    ### define plot

    fig = go.Figure(
        data=[
            entropy_noise_binary_trace,
            entropy_noise_star_trace,
            phylogenetic_noise_trace
        ],
        layout=go.Layout(
            title="quantification of noise",
            font=dict(size=18),
            margin=dict(b=45, t=50),
            xaxis=dict(
                title="#predicted contacts / protein length",
                showspikes=True
            ),
            yaxis=dict(
                title="fraction of noise",
                range=[0,0.8],
                showspikes=True
            ),
            legend=dict(
            orientation="v",
            x=0.15, y=1.0
            ),
            width=width,
            height=height
        )
    )

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


    args = parse_args()

    data_dir = args.data_dir
    plot_dir = data_dir+"/plots/benchmarks/"
    pdb_dir = data_dir+"/pdb/"

    sequence_separation = 6
    contact_thr = 8
    non_contact_thr = 8


    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)



    ### create benchmark plot for star-tree topologies

    # create benchmark object
    b = Benchmark(pdb_dir)

    #specify methods to benchmark
    b.add_method("APC", data_dir +"/recover_pcd_constrained/", "apc.star.mat")
    b.add_method("EC", data_dir +"/recover_pcd_constrained/", "ec.star.mat")
    b.add_method("no APC", data_dir + "/recover_pcd_constrained/", "raw.star.mat")

    #add constraint that all MRF optimizations have exist status 0
    b.add_constraint("opt_code", 0, "greater_equal")

    #compute the precision of predictions
    b.compute_evaluation_statistics(seqsep=sequence_separation, contact_thr=contact_thr, noncontact_thr=non_contact_thr)

    #generate a benchmark plot
    benchmark_plot_star = b.plot_precision_vs_rank()
    plot_file = plot_dir+"/"+"fig_6b.html"
    plot_ccmgen_benchmark_figure(benchmark_plot_star, 'star topology', plot_file, height=350, width=500)






    ### create benchmark plot for binary-tree topologies

    # create benchmark object
    b = Benchmark(pdb_dir)

    # specify methods to benchmark
    b.add_method("APC", data_dir + "/recover_pcd_constrained/", "apc.binary.mat")
    b.add_method("no APC", data_dir + "/recover_pcd_constrained/", "raw.binary.mat")
    b.add_method("EC", data_dir + "/recover_pcd_constrained/", "ec.binary.mat")

    # add constraint that all MRF optimizations have exist status 0
    b.add_constraint("opt_code", 0, "greater_equal")

    # compute the precision of predictions
    b.compute_evaluation_statistics(seqsep=sequence_separation, contact_thr=contact_thr, noncontact_thr=non_contact_thr)

    # generate a benchmark plot
    benchmark_plot_binary = b.plot_precision_vs_rank()
    plot_file = plot_dir+"/"+"fig_6a.html"
    plot_ccmgen_benchmark_figure(benchmark_plot_binary, 'binary topology', plot_file, height=350, width=500)






    ### create qunatification of noise plot
    plot_file = plot_dir+"/"+"fig_6c.html"
    plot_ccmgen_noise_quant_figure(benchmark_plot_star, benchmark_plot_binary, plot_file, height=350, width=500)





if __name__ == '__main__':
    main()