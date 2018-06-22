#!/usr/bin/env bash

#------------------------------------------------------------------------------
# This script will reproduce Figures in style of 3d
#   it requires alignments from the PSICOV dataset and
#   binary raw files for MRF models generated in step 1 and 2
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------

data_dir=$1

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
plot_dir=$data_dir"/plots/apc_vs_ec/"
alignment_dir=$data_dir"/aln/"
predictions_pll=$data_dir"/predictions_pll/"
predictions_pcd=$data_dir"/predictions_pcd/"

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
    echo "Generating scatter plots for the APC term vs the entropy correction term for protein $name..."


    #generate plot for MRF model learned by maximizing pseudo-likelihood
    binary_raw_pll=$predictions_pll"/$name.braw.gz"
    if [ -f $binary_raw_pll ]
    then
        plot_name_pll=$plot_dir"/$name.apc_vs_ec.pll.html"
        python $script_dir/../plotting/scripts/plot_scatter_apc_vs_ec.py $binary_raw_pll $alignment_file $plot_name_pll
    fi


    #generate plot for MRF model learned with persistent constrastive divergence
    binary_raw_pcd=$predictions_pcd"/$name.braw.gz"
    if [ -f $binary_raw_pcd ]
    then
        plot_name_pcd=$plot_dir"/$name.apc_vs_ec.pcd.html"
        python $script_dir/../plotting/scripts/plot_scatter_apc_vs_ec.py $binary_raw_pcd $alignment_file $plot_name_pcd
    fi

done