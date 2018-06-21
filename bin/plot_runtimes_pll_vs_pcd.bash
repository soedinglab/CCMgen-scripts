#!/usr/bin/env bash

#------------------------------------------------------------------------------
# This script will reproduce Figure 1d
#   requires contact matrix files generated in step 1 and 2
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------

data_dir=$1

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
plot_dir=$data_dir"/plots/"
pll_dir=$data_dir"/predictions_pll/"
pcd_dir=$data_dir"/predictions_pcd/"

if [ ! -d $plot_dir ]
then
    mkdir $plot_dir
fi

#------------------------------------------------------------------------------
# run the plotting script
#------------------------------------------------------------------------------

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python $script_dir/../plotting/scripts/plot_runtimes_pll_vs_pcd.py \
    $plot_dir \
    "$pll_dir,$pcd_dir" \
    "pseudo-likelihood,persistent contrastive divergence"