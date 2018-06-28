#!/usr/bin/env bash

#------------------------------------------------------------------------------
# This script will reproduce the contact map plots in Figure 3
#   It requires matrix files that have been generated in step 1 and 2
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------

data_dir=$1

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
plot_dir=$data_dir"/plots/contact_maps/"
mat_dir_pll=$data_dir"/predictions_pll/"
mat_dir_pcd=$data_dir"/predictions_pcd/"
alignment_dir=$data_dir"/aln/"


if [ ! -d $plot_dir ]
then
    mkdir -p $plot_dir
fi


#------------------------------------------------------------------------------
# run the plotting script for pseudo-likelihood
#------------------------------------------------------------------------------

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for alignment_file in $(ls $alignment_dir/*.aln);
do
    name=$(basename $alignment_file .aln)
    echo "Plotting contact maps for protein $name..."


    #generate plot for MRF model learned by maximizing pseudo-likelihood
    raw_mat_pll=$mat_dir_pll"/$name.mat"

    if [ -f $raw_mat_pll ]
    then
        plot_name=$plot_dir"/$name.contact_map.pll.html"
        ccm_plot cmap --mat-file $raw_mat_pll  -o $plot_name --seq-sep 1 --contact-threshold 8

        plot_name=$plot_dir"/$name.contact_map.pll.apc.html"
        ccm_plot cmap --mat-file $raw_mat_pll  -o $plot_name --seq-sep 1 --contact-threshold 8 \
        --apc

        raw_mat_ec_pll=$mat_dir_pll"/$name.ec.mat"
        if [ -f $raw_mat_ec_pll ]
        then
            plot_name=$plot_dir"/$name.contact_map.pll.ec.html"
            ccm_plot cmap --mat-file $raw_mat_ec_pll  -o $plot_name  --seq-sep 1 --contact-threshold 8
        fi
    fi

    #generate plot for MRF model learned by maximizing pseudo-likelihood
    raw_mat_pcd=$mat_dir_pcd"/$name.mat"

    if [ -f $raw_mat_pcd ]
    then
        plot_name=$plot_dir"/$name.contact_map.pcd.html"
        ccm_plot cmap --mat-file $raw_mat_pcd  -o $plot_name --seq-sep 1 --contact-threshold 8

        plot_name=$plot_dir"/$name.contact_map.pcd.apc.html"
        ccm_plot cmap --mat-file $raw_mat_pcd  -o $plot_name --seq-sep 1 --contact-threshold 8 --apc

        raw_mat_ec_pcd=$mat_dir_pcd"/$name.ec.mat"
        if [ -f $raw_mat_ec_pcd ]
        then
            plot_name=$plot_dir"/$name.contact_map.pcd.ec.html"
            ccm_plot cmap --mat-file $raw_mat_ec_pcd  -o $plot_name  --seq-sep 1 --contact-threshold 8
        fi
    fi


done