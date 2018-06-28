#!/usr/bin/env bash


#------------------------------------------------------------------------------
# This script runs CCMpredPy with pseudo-likelihood maximization for all
#   Pfam alignments in the PSICOV dataset.
#   The path to the directory containing the PSICOV data needs to be specified
#   by the first argument.
#   The second arguments specifies the number of OMP threads for parallelization.
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# parameters
#------------------------------------------------------------------------------
data_dir=$1
num_threads=$2

#------------------------------------------------------------------------------
# set up OpenMP
#------------------------------------------------------------------------------
export OMP_NUM_THREADS=$num_threads
echo "using " $OMP_NUM_THREADS "threads for omp parallelization"

#------------------------------------------------------------------------------
# create data structure
#------------------------------------------------------------------------------
alignment_dir=$data_dir"/aln/"
mat_dir=$data_dir"/predictions_pll/"

if [ ! -d $mat_dir ]
then
    mkdir $mat_dir
fi


#------------------------------------------------------------------------------
# settings for CCMpredPy
#------------------------------------------------------------------------------
settings=" --num-threads $num_threads --aln-format psicov"
settings=$settings" --wt-simple"
settings=$settings" --max-gap-seq 75 --max-gap-pos 50"
settings=$settings" --reg-lambda-single 10 --reg-lambda-pair-factor  0.2 --v-center"
settings=$settings" --pc-uniform --pc-single-count 1 --pc-pair-count 1"
settings=$settings" --ofn-pll "
settings=$settings" --maxit 5000 --lbfgs-ftol 1e-4 --lbfgs-max-linesearch 5 --lbfgs-maxcor 5"


#------------------------------------------------------------------------------
# Run CCMpredPy
#------------------------------------------------------------------------------

for alignment_file in $(ls $alignment_dir/*.aln);
do

    name=$(basename $alignment_file ".aln")

    file_paths=" -b "$mat_dir"/"$name".braw.gz "
    file_paths=$file_paths" -m $mat_dir/$name.raw.mat "
    file_paths=$file_paths" --apc $mat_dir/$name.apc.mat "
    file_paths=$file_paths" --entropy-correction $mat_dir/$name.ec.mat "
    file_paths=$file_paths" $alignment_file "
    log_file=$mat_dir"/"$name.log

    echo -e "Running CCMpredPy with pll maximization for protein: $name \n(Status is logged in: $log_file)"

    if [ ! -f $log_file ]
    then
        ccmpred $settings" "$file_paths  > $log_file
    fi
done