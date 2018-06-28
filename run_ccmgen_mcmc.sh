#!/usr/bin/env bash


#------------------------------------------------------------------------------
# This script runs CCMgen  to create MCMC samples from a pre-defined Markov
# random field model that is specified by a .braw.gz file.
#   The path to the directory containing the binary raw files that define the
#   Markov random field model needs to be specified with the first argument.
#   The second arguments specifies the number of OMP threads for parallelization.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------
binary_raw_dir=$1
num_threads=$2

#------------------------------------------------------------------------------
# set up OpenMP
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=$num_threads
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
algorithm=$(echo $(basename $binary_raw_dir) | sed -e "s/^predictions_//")
sample_dir=$(dirname $binary_raw_dir)"/samples_$algorithm"
alignment_dir=$(dirname $binary_raw_dir)"/aln/"

if [ ! -d $sample_dir ]
then
    mkdir $sample_dir
fi

#------------------------------------------------------------------------------
# settings for CCMgen
#------------------------------------------------------------------------------

settings=" --aln-format psicov --max-gap-pos 50 --max-gap-seq 75 --num-threads $num_threads"
settings=$settings" --mcmc-sampling --mcmc-sample-random-gapped --mcmc-burn-in 500 --num-sequences 10000"


#------------------------------------------------------------------------------
# run CCMgen
#------------------------------------------------------------------------------

for alignment_file in $(ls $alignment_dir/*.aln);
do

    name=$(basename $alignment_file .aln)

    braw_file="$binary_raw_dir/$name.braw.gz"
    sample_file="$sample_dir/$name.mcmc.aln"
    log_file="$sample_dir/$name.mcmc.log"

    if [ ! -f $log_file ] && [ -f $braw_file ]
    then

        file_paths=" --alnfile "$alignment_file
        file_paths=$file_paths" "$braw_file" "$sample_file

        echo -e "Running CCMgen for MRF learned with $algorithm to generate MCMC sample for protein: $name \n(Status is logged in: $log_file)"
        ccmgen $settings" "$file_paths > $log_file
    fi
done