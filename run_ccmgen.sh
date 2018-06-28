#!/usr/bin/env bash


#------------------------------------------------------------------------------
# This script runs CCMgen to create a synthetic alignment that has been sampled
# along a STAR tree from a pre-defined Markov random field model that is
# specified by a .braw.gz file.
# The path to the directory containing the binary raw files that define the
# Markov random field model needs to be specified with the first argument.
# The second arguments specifies the number of OMP threads for parallelization.
# The third argument specifies the topology along which new samples are
# generated (either star or binary).
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------
data_dir=$1
num_threads=$2
topology=$3

#------------------------------------------------------------------------------
# set up OpenMP
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=$num_threads
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
binary_raw_dir=$data_dir"/predictions_pcd_constrained/"
sample_dir=$data_dir"/samples_pcd_constrained"
alignment_dir=$data_dir"/aln/"

if [ ! -d $sample_dir ]
then
    mkdir $sample_dir
fi

#------------------------------------------------------------------------------
# settings for CCMgen
#------------------------------------------------------------------------------

settings=" --aln-format psicov --max-gap-pos 50 --max-gap-seq 75 --num-threads $num_threads"
settings=$settings" --tree-$topology --mutation-rate-neff --burn-in 10"

#------------------------------------------------------------------------------
# run CCMgen
#------------------------------------------------------------------------------

for alignment_file in $(ls $alignment_dir/*.aln);
do

    name=$(basename $alignment_file .aln)

    braw_file="$binary_raw_dir/$name.braw.gz"
    sample_file="$sample_dir/$name.$topology.aln"
    log_file="$sample_dir/$name.$topology.log"

    if [ ! -f $log_file ] && [ -f $braw_file ]
    then

        file_paths=" --alnfile "$alignment_file
        file_paths=$file_paths" "$braw_file" "$sample_file

        echo -e "Running CCMgen for MRF learned with PCD (constrained) to generate
        \n\t synthetic alignment ($topology topology) for protein: $name
        \n\t (Status is logged in: $log_file)"

        ccmgen $settings" "$file_paths > $log_file
    fi
done