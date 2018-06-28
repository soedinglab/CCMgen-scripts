from collections import Counter
import glob
import os
import ccmpred.io.pdb
import ccmpred.io.contactmatrix
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
import colorlover as cl


class Benchmark():
    """
    Benchmarking contact prediction methods on a dataset
    """

    def __init__(self, pdb_dir):
        self.pdb_dir = pdb_dir
        self.pdb_files = glob.glob(self.pdb_dir +"/*pdb")

        self.evaluation_data = {}
        self.ordered_methods = []
        self.evaluation_statistics = {}
        self.filter = []

    def __apply_filter(self, protein):

        filter_operators = {
            'greater': np.greater,
            'less': np.less,
            'greater_equal': np.greater_equal,
            'less_equal': np.less_equal,
            'equal': np.equal,
            'not_equal': np.not_equal
        }

        for method_name in self.ordered_methods:

            mat_file = self.evaluation_data[method_name][protein]
            mat, meta = ccmpred.io.contactmatrix.read_matrix(mat_file)

            for f in self.filter:
                filter_res = ccmpred.io.contactmatrix.find_dict_key(f['key'], meta)
                if not filter_operators[f['operator']](filter_res, f['value']):
                    print("{0} did not pass filter for {1} {2} {3}: {4}".format(method, f['key'], f['operator'], f['value'], filter_res))
                    return False

        return True

    def __get_distances(self, pdb_file, L):

        if not os.path.exists(pdb_file):
            raise IOError("PDB File " + str(pdb_file) + "does not exist. ")

        # determine distance matrix from PDB file
        distance_matrix = ccmpred.io.pdb.distance_map(pdb_file, L=L)

        # get residue pairs that are resolved (not NAN)
        index_i, index_j = np.where(~np.isnan(distance_matrix))

        # Create the evaluation file
        eval_df = pd.DataFrame(
            {
                'i': index_i,
                'j': index_j,
                'cb_distance': distance_matrix[index_i, index_j],
            }
        )
        eval_df.sort_values(by=['i', 'j'], inplace=True)

        #do not count residue pairs twice
        eval_df = eval_df[eval_df['j'] > eval_df['i']]

        return eval_df

    def __compute_precision_recall(self, true_class, score):
        """
        Compute Precision and Recall for a running threshold

        :param true_class: binary vector
        :param score: score that is used for ranking
        :return:
        """

        df = pd.DataFrame({'true':true_class, 'score':score})
        df.sort_values('score', ascending=False, inplace=True)

        df['cumsum_pred']   = range(1, len(df)+1)
        df['cumsum_tp']     = df.true.cumsum()


        df['precision'] = df['cumsum_tp']  / df['cumsum_pred']
        df['recall']    = df['cumsum_tp']  / np.sum(true_class)

        return df.precision.tolist(), df.recall.tolist(), df.score.tolist()

    def __compute_mean_error(self,cb_distance, score, contact_thr):
        """
            Compute mean error of predictions:

            error:  0, iff cb_distance <= contact_thr
            d, d=cb_distance - contact_thr

        :param cb_distance:
        :param score:
        :param contact_thr:
        :return:
        """

        df = pd.DataFrame({'cb_distance':cb_distance, 'score':score})
        df.sort_values('score', ascending=False, inplace=True)

        #compute the error
        df.loc[:,'error'] = df.cb_distance - contact_thr
        df.loc[df.error < 0 ,'error'] = 0

        df.loc[:,'mean_error'] = df.error.expanding(1).mean()

        return df.mean_error.tolist()

    def __compute_evaluation_statistics_protein(self, pdb_file, ranks, seqsep, contact_thr, noncontact_thr, meta):

        protein = os.path.basename(pdb_file).split(".")[0]

        # compute distance matrix from PDB file
        eval_df = self.__get_distances(pdb_file, meta['L'])

        # remove pairs that are separated less than SEQSEP positions along primary sequence
        eval_df['class'] = (eval_df['cb_distance'] <= contact_thr) * 1
        eval_df = eval_df[eval_df['j'] >= (eval_df['i'] + seqsep)]

        # in case noncontact_thr != contact_thr: remove residue pairs with contact_thr < Cb distance < noncontact_thr
        if noncontact_thr > contact_thr:
            eval_df = eval_df[(eval_df['cb_distance'] <= contact_thr) | (eval_df['cb_distance'] > noncontact_thr)]

        #sort residue pairs
        eval_df.sort_values(by=['i', 'j'], inplace=True)
        eval_df.reset_index(inplace=True, drop=True)

        # add scores from all methods
        for method_name in self.ordered_methods:
            mat, _ = ccmpred.io.contactmatrix.read_matrix(self.evaluation_data[method_name][protein])
            eval_df[method_name] = mat[eval_df['i'], eval_df['j']]

        # determine number of top ranked residue pairs that will be considered for evaluation
        ranks_L = np.round(meta['L'] * ranks).astype(int)
        # if there are less residue pairs than max(rank_L): adjust rank_L
        ranks_L = np.array([rank for rank in ranks_L if rank < len(eval_df)])


        # compute precision and recall values for all methods
        protein_eval_metrics={}
        for method_name in self.ordered_methods:

            precision, recall, threshold = self.__compute_precision_recall(
                eval_df['class'], eval_df[method_name])
            mean_error = self.__compute_mean_error(eval_df['cb_distance'], eval_df[method_name], contact_thr)

            protein_eval_metrics[method_name] = {}
            protein_eval_metrics[method_name]['precision'] = [np.array(precision)[rank] for rank in ranks_L]
            protein_eval_metrics[method_name]['mean_error'] = [np.array(mean_error)[rank] for rank in ranks_L]
            protein_eval_metrics[method_name]['recall'] = [np.array(recall)[rank] for rank in ranks_L]

        return protein_eval_metrics

    def __compute_meanprecision_per_rank(self):

        mean_precision_per_rank = {}
        mean_precision_per_rank['ranks'] = self.evaluation_statistics['ranks']
        for method_name in self.ordered_methods:

            precision_per_protein = [self.evaluation_statistics['proteins'][protein][method_name]['precision'] for protein in self.evaluation_statistics['proteins']]
            mean_precision_per_rank[method_name] = np.nanmean(precision_per_protein, axis=0)

        return mean_precision_per_rank

    def __plot_precision_vs_rank_plotly(self, mean_precision_per_rank, title, yaxistitle, legend_order=None, plot_file=None):

        # define order of methods in the legend
        methods = legend_order
        if legend_order is None:
            methods = list(mean_precision_per_rank.keys())
            methods.remove('ranks')

        # define suitable colors
        method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set1'])
        if 10 < len(methods) < 13:
            method_colors = np.array(cl.scales[str(max(3, len(methods)))]['qual']['Set3'])
        elif len(methods) >= 13:
            method_colors = np.array(cl.to_rgb(cl.interp(cl.scales['9']['qual']['Set1'], 50)))


        data = []
        for nr, method in enumerate(methods):

            method_trace = go.Scatter(
                x=[str(rank) for rank in np.round(mean_precision_per_rank['ranks'], decimals=2)],
                y=mean_precision_per_rank[method],
                name=method,
                mode='lines',
                line=dict(
                    width=4,
                    color=method_colors[nr]
                )
            )
            data.append(method_trace)


        plot = {
            "data": data,
            "layout": go.Layout(
                hovermode = 'closest',
                title=title,
                xaxis1=dict(
                    title='#predicted contacts / protein length',
                    tickvals=[str(rank) for rank in np.linspace(1, 0, 10, endpoint=False)[::-1]]),
                yaxis1=dict(
                    title=yaxistitle,
                    range=[0, 1]
                ),
                font=dict(size=18)
            )
        }

        #move plot a bit upwards if there is no title
        if title=="":
            plot["layout"]['margin']['t'] = 10

        if plot_file is not None:
            plotly_plot(plot, filename=plot_file, auto_open=False, show_link=False)
        else:
            return plot

    def add_method(self, method_name, method_dir, filter=""):

        self.ordered_methods.append(method_name)

        self.evaluation_data[method_name] = {}

        for mat_file in glob.glob(method_dir + "/*" + filter + "*"):
            protein = os.path.basename(mat_file).split(".")[0]
            self.evaluation_data[method_name][protein] = mat_file

    def reset_methods(self):
        self.ordered_methods = []
        self.evaluation_data = {}
        self.evaluation_statistics = {}
        self.filter = []

    def add_constraint(self, key, value, operator):
        self.filter.append(
            {"key": key,
             "value": value,
             "operator": operator}
        )

    def compute_evaluation_statistics(self, seqsep=12, contact_thr=8, noncontact_thr=8):

        self.evaluation_statistics = {}

        # definition of true positive (residue-residue contact based on distance between Cb atoms)
        self.evaluation_statistics['contact_thr'] = contact_thr

        # definition of true negative (what is NOT a residue-residue contact based on distance between Cb atoms)
        if noncontact_thr < contact_thr:
            noncontact_thr = contact_thr
        self.evaluation_statistics['noncontact_thr'] = noncontact_thr

        # ignore residue pairs that are separated by less than SEQSEP positions in the primary sequence
        self.evaluation_statistics['seqsep'] = seqsep

        # define x-axis: number of top ranked predictions (wrt to protein length) that will be considered for evaluation
        self.evaluation_statistics['ranks'] =  np.linspace(1, 0, 50, endpoint=False)[::-1]

        # name of methods for evaluation
        self.evaluation_statistics['methods'] = self.ordered_methods

        #dictionary collecting the evaluation statistics per protein
        self.evaluation_statistics['proteins'] = {}

        print("Compute evaluation statistics for {0} proteins and methods:".format(len(self.pdb_files)))
        print(self.ordered_methods)

        # iterate over proteins with pdb structures
        for id, pdb_file in enumerate(self.pdb_files):

            protein = os.path.basename(pdb_file).split(".")[0]
            print(str(id + 1) + "/" + str(len(self.pdb_files)) + " " + str(protein))

            # ensure that all methods are compared on the same data set:
            # if a protein is not available for one of the methods then this protein is skipped
            if not all([protein in self.evaluation_data[method_name].keys() for method_name in
                        self.ordered_methods]):
                print("No scores available for protein {0} for at least one of the methods".format(protein))
                continue

            # ensure that special constraints are fulfilled
            if not (self.__apply_filter(protein)):
                print("Protein {0} did not pass filters for at least one of the methods.".format(protein))
                continue

            #get some meta information about the protein from one of the methods meta info
            meta_protein = {}
            mat_file = self.evaluation_data[list(self.evaluation_data.keys())[0]][protein]
            mat, meta = ccmpred.io.contactmatrix.read_matrix(mat_file)
            meta_protein['L'] = ccmpred.io.contactmatrix.find_dict_key('ncol', meta)
            meta_protein['N'] = ccmpred.io.contactmatrix.find_dict_key('nrow', meta)
            meta_protein['Diversity'] = ccmpred.io.contactmatrix.find_dict_key('diversity', meta)
            meta_protein['neff'] = ccmpred.io.contactmatrix.find_dict_key('neff', meta)

            # compute evaluation metrics: precision, recall, mean error for every method in benchmark_methods
            self.evaluation_statistics['proteins'][protein] = self.__compute_evaluation_statistics_protein(
                pdb_file, self.evaluation_statistics['ranks'], seqsep, contact_thr, noncontact_thr, meta_protein)

        print("There are {0} proteins in the evaluation data set.".format(len(self.evaluation_statistics['proteins'])))

    def plot_precision_vs_rank(self, plot_file=None):

        if len(self.evaluation_statistics) == 0:
            print("You first need to calculate statistics for selected methods!")
            return

        mean_precision_per_rank = self.__compute_meanprecision_per_rank()

        title=""
        yaxistitle = 'Mean Precision over Proteins'

        return self.__plot_precision_vs_rank_plotly(
                mean_precision_per_rank, title, yaxistitle, legend_order=self.ordered_methods, plot_file=plot_file)