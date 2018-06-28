# CCMgen Scripts

A collection of scripts to reproduce the results presented in the publication: 
"Synthetic protein alignments by CCMgen quantify noise in residue-residue contact prediction", currently available at bioRxiv, doi: [10.1101/344333](https://doi.org/10.1101/344333).


## Installation

Clone this repository using:

```bash
git clone https://github.com/soedinglab/ccmgen-scripts.git
```

Then, change into the directory and install the python packages required for plotting by running the command:

```bash
pip install .
```

## Download Data Set

First of all, download the [PSICOV dataset](http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/suppdata/) which is the basis for the following analysis.


## Reproduce Raw Data

Next, specify the path to the PSICOV data set that you just downloaded and the number of OMP threads that you want to use for parallelization:

```bash
data_dir=path/to/psicov/data/
num_threads=4

```

In order to generate all the necessary data to reproduce the analysis presented in the paper, run the following scripts from the CCMgen-scripts repository main folder:

1. ```bash run_ccmpred_pll.sh $data_dir $num_threads```

	This command will call CCMpredPy and predict contact maps by maximizing the pseudo-likelihood. 
	It will generate *.mat, *.apc.mat and *.braw.gz files in a folder named "predictions_pll" for all proteins in the data set that are required to reproduce Figure 1.

2. ```bash run_ccmpred_pcd.sh $data_dir $num_threads```

      This command will call CCMpredPy and predict contact maps by learning a Markov random field with persistent contrastive divergence. 
      It will generate *.mat, *.apc.mat and *.braw.gz files in a folder named "predictions_pcd" for all proteins in the data set that are required to reproduce Figure 1.

3. ```bash run_ccmgen_mcmc.sh $data_dir/predictions_pll $num_threads```

	This command will call CCMgen and generate 10000 protein sequences from a pre-defined Markov random field (learned in step 1) via MCMC sampling.
	It will generate *.aln files in a new directory named samples_pll.

4. ```bash run_ccmgen_mcmc.sh $data_dir/predictions_pcd $num_threads```

      This command will call CCMgen and generate 10000 protein sequences from a pre-defined Markov random field (learned in step 2) via MCMC sampling.
      It will generate *mcmc.aln files in a new directory named samples_pcd.

5. ```bash run_ccmpred_pcd_with_constraints.sh $data_dir $num_threads```

	This command will call CCMpredPy and predict contact maps by learning a Markov random field with persistent contrastive divergence.
	Residue pairs that do not form contacts (C_beta distance > 12 angstrom) in the reference protein structure will receive zero couplings. 
	Running this command will generate *.mat and *.braw.gz files in a folder named "predictions_pcd_constrained" for all proteins in the data set.

6. ```bash run_ccmgen.sh $data_dir $num_threads star```

	This command will call CCMgen to generate a synthetic alignment along a STAR-tree topology and according to the constraints from Markov random fields that have been learned in step 5. 
	It will generate *star.aln files in a new directory named "samples_pcd_constrained" for all proteins in the data set.

7. ```bash run_ccmgen.sh $data_dir $num_threads binary```
        This command will call CCMgen to generate a synthetic alignment along a BINARY-tree topology and according to the constraints from Markov random
        fields that have been learned in step 5.
        It will generate *binary.aln files in a new directory named "samples_pcd_constrained" for all proteins in the data set.
	
8. ```bash run_ccmpred_pcd_recover.sh $data_dir $num_threads star```

	This command will call CCMpredPy and predict contact maps by learning a Markov random field with persistent contrastive divergence.
	The input alignments are the synthetic alignments generated with CCMgen along a star-tree topology in step 6.
	It will generate *.mat, *.apc.star.mat and *.ec.star.mat files in a folder named "recover_pcd_constrained" for all proteins in the data set.


9. ```bash run_ccmpred_pcd_recover.sh $data_dir $num_threads binary```

      This command will call CCMpredPy and predict contact maps by learning a Markov random field with persistent contrastive divergence.
      The input alignments are the synthetic alignments generated with CCMgen along a binary-tree topology in step 7.
      It will generate *.mat, *.apc.binary.mat and *.ec.binary.mat files in a folder named "recover_pcd_constrained" for all proteins in the data set.

## Reproduce Figure 1

1. ```bash plot_fig_1ab.sh $data_dir```

	This command will generate plots like Figure 1A and 1A for all proteins in the PSICOV data set.
	In order to generate the plots, MRF models need to be learned by maximizing pseudo-likelihood and persistent contrastive divergence as described in step 1 and 2. 
	Furthermore, MCMC samples need to be generated in advance as described in step 3 and 4.
	Note: the generated .html files can become large!
	Plots will be written to ```$data_dir/plots/alignment_statistics/```.

2. ```python plot_fig_1c.py $data_dir```

	This script will reproduce the contact prediction benchmark for pseudo-likelihood and persistent contrastive divergence from Figure 1C. 	
	In order to generate the plots, MRF models need to be learned by maximizing pseudo-likelihood and persistent contrastive divergence as described in step 1 and 2.
	The plot is written to ```$data_dir/plots/benchmarks/fig_1c.html```.
	
3. ```python plot_fig_1d.py $data_dir```

	This command will generate boxplots to visualize the distribution of run times for learning MRF with CCMpredPy like in Figure 1D for all proteins in the data set. 
	In order to generate this plot, MRF models need to be learned by maximizing pseudo-likelihood and persistent contrastive divergence as described in step 1 and 2.
	The plot will be written to ```$data_dir/plots/fig_1d.html```.

## Reproduce Figure 3

1. ```bash plot_fig_3abc.sh $data_dir```

	This command will generate plots of the contact score matrices, comprising the raw contact scores and APC and entropy corrected scores for all proteins in the PSICOV data set.
	In order to generate the plots, MRF models need to be learned by maximizing pseudo-likelihood and persistent contrastive divergence as described in step 1 and 2.
	Plots will be written to ```$data_dir/plots/contact_maps/```.

2. ```python plot_fig_3d.py $data_dir```

	This command will generate scatter plots of the APC correction term vs the entropy correction term per residue pair for all proteins in the PSICOV data set.
	In order to generate the plots, MRF models need to be learned by maximizing pseudo-likelihood and persistent contrastive divergence as described in step 1 and 2.
	Plots will be written to ```$data_dir/plots/apc_vs_ec/```.

## Reproduce Figure 6



