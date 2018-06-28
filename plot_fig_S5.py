#!/usr/bin/env python

# ===============================================================================
###     This script plots Figure Supplemental Figure 5 of the CCMgen paper:
###     boxplots visualizing the distribution of mutation rates used for
###     generating the synthetic alignments with CCMgen
###     and boxplots visualizing the difference in Neff values between
###     synthetic and original Pfam alignments
# ===============================================================================

### load libraries
import argparse
import os
from benchmark import Benchmark
import copy
from plotly.offline import plot as plotly_plot
import plotly.graph_objs as go
import glob

def plot_boxplot(statistics_dict, property, plot_file):

    topologies = sorted(statistics_dict.keys())

    data = []
    for topology in topologies:

        values = statistics_dict[topology][property]
        proteins = statistics_dict[topology]['protein']
        target_neff = statistics_dict[topology]['target neff']
        sample_neff = statistics_dict[topology]['sample neff']

        hover_text = ["{0}<br>target neff:{1}<br>sample neff:{2}".format(
            proteins[i], target_neff[i], sample_neff[i]) for i in range(len(values))]

        box = go.Box(
            y=values,
            boxmean=True,
            pointpos=1.8,
            jitter=0.4,
            boxpoints='all',
            name=topology,
            marker=dict(opacity=1),
            text=hover_text,
            hoverinfo='all',
            orientation='v',
            showlegend=False
        )

        data.append(box)


    plot = {
        "data": data,
        "layout": go.Layout(
            yaxis=dict(
                exponentformat='e',
                showexponent='All'
            ),
            font=dict(size=18)
        )
    }

    if property == "neff_difference":
        plot['layout']['yaxis']['title'] = "Pfam Neff - synthetic Neff"
    if property == "mutation_rate":
        plot['layout']['yaxis']['title'] = "mutation rate"

    plotly_plot(plot, filename=plot_file, auto_open=False)

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
    sampled_aln = data_dir+"/samples_pcd_constrained/"

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)


    statistics_dict={
        'star':{
            'protein' : [],
            'neff_difference' : [],
            'target neff': [],
            'sample neff': [],
            'mutation_rate': []
        },
        'binary':{
            'protein': [],
            'neff_difference': [],
            'target neff': [],
            'sample neff': [],
            'mutation_rate': []
        }
    }

    log_files = glob.glob(sampled_aln + "/*.log")
    for log_file in log_files:

        topology = os.path.basename(log_file).split(".")[-2]
        protein = os.path.basename(log_file).split(".")[0]

        #read log file
        with open(log_file) as f:
            content = f.readlines()

        if len(content) == 0:
            print("no content", log_file)
            continue

        #parse original Pfam Neff
        target_neff_list = [line for line in content if "Neff(HHsuite-like)=" in line]
        if len(target_neff_list) == 0:
            print("no Neff", log_file)
            continue
        target_neff = float(target_neff_list[0].split("=")[2].replace(".\n", ""))


        # parse latest sampled Neff
        sampled_neff_list = [line for line in content if "has Neff" in line]
        if len(sampled_neff_list) == 0 :
            print("no sample Neff", log_file)
            continue

        sampled_neff = float(sampled_neff_list[-1].split(" ")[13])
        diff=target_neff - sampled_neff


        #parse latest mutation rate
        mutation_rate_list = [line for line in content if "mutation rate" in line]
        mutation_rate = mutation_rate_list[-1].split(" ")[9]


        statistics_dict[topology]['protein'].append(protein)
        statistics_dict[topology]['neff_difference'].append(diff)
        statistics_dict[topology]['target neff'].append(target_neff)
        statistics_dict[topology]['sample neff'].append(sampled_neff)
        statistics_dict[topology]['mutation_rate'].append(mutation_rate)



    plot_boxplot(statistics_dict, "neff_difference", plot_dir + "/fig_S5a.html")
    plot_boxplot(statistics_dict, "mutation_rate", plot_dir + "/fig_S5b.html")





if __name__ == '__main__':
    main()