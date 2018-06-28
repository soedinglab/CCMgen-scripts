#!/usr/bin/env bash

#------------------------------------------------------------------------------
# This script will reproduce Figures in style of 1a and 1b
#   requires alignments from PSICOV dataset and
#   MCMC samples generated from MRF models generated in step 1,2 and 3
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------

data_dir=$1

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
plot_dir=$data_dir"/plots/alignment_statistics/"
alignment_dir=$data_dir"/aln/"
samples_pll_dir=$data_dir"/samples_pll/"
samples_pcd_dir=$data_dir"/samples_pcd/"

if [ ! -d $plot_dir ]
then
    mkdir -p $plot_dir
fi


#------------------------------------------------------------------------------
# run the plotting script
#------------------------------------------------------------------------------

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for alignment_file in $(ls $alignment_dir/*.aln);
do
    name=$(basename $alignment_file .aln)
    echo "Plotting alignment statistiscs for protein $name..."


    #generate plot for MRF model learned by maximizing pseudo-likelihood
    sampled_alignment_pll=$samples_pll_dir"/$name.mcmc.aln"

    if [ -f $sampled_alignment_pll ]
    then
        plot_name_pll=$plot_dir"/$name.alignment_stats_mcmc_vs_observed.pll.html"
        ccm_plot cmap -a $alignment_file -s $sampled_alignment_pll -o $plot_name_pll
    fi


    #generate plot for MRF model learned with persistent constrastive divergence
    sampled_alignment_pcd=$samples_pcd_dir"/$name.mcmc.aln"
    if  [ -f $sampled_alignment_pcd ]
    then
        plot_name_pcd=$plot_dir"/$name.alignment_stats_mcmc_vs_observed.pcd.html"
        ccm_plot aln-stats -a $alignment_file -s $sampled_alignment_pcd -o $plot_name_pcd
    fi

done