import sys
sys.path.insert(0, "..")

import os

import numpy as np
import pandas as pd

from scipy.stats import beta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['figure.figsize'] = 15, 10

from config.config import config, read_config
from config.parse_input_args import parse_args
from special_utils.os_utils import smart_makedirs

from chains.chains import Chains
from kmer_utilities.kmer_utilities import kmer_names
from mutation_strategies.mutation_strategies import MutationStrategies

from model_analysis import read_models

import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot


def plot_distr(models, figures_dir):
    def plot_distr_str(models, fig_dir, beta_step=0.001):
        betas_fr, betas_cdr, betas_full, dirs_subst = {}, {}, {}, {}
        for chain, model in models.iteritems():
            betas_fr[chain] = np.array(model[["beta_FR_shape1",
                                              "beta_FR_shape2"]])
            betas_cdr[chain] = np.array(model[["beta_CDR_shape1",
                                               "beta_CDR_shape2"]])
            betas_full[chain] = np.array(model[["beta_FULL_shape1",
                                                "beta_FULL_shape2"]])
            dirs_subst[chain] = np.array(model[["dir_shape1",
                                                "dir_shape2",
                                                "dir_shape3"]])
        def add_distr(params, chain, i, segment):
            param = params[chain][i, :]
            if np.any(np.isnan(param)):
                return None
            support = np.arange(0, 1, beta_step)
            pdf = beta.pdf(support, *param)
            trace = go.Scatter(
                x = support,
                y = pdf,
                mode = 'lines',
                name = chain.name + " " + segment
            )
            return trace

        def append_notnan(l, o):
            if o is not None:
                l.append(o)

        def plot_fig(data, layout, filename):
            if data:
                fig = dict(data=data, layout=layout)
                plot(fig, filename=filename, auto_open=False)

        for i, kmer in enumerate(kmer_names()):
            print(kmer)
            grs_fr, grs_cdr, grs_full = [], [], []
            layout = dict(title = 'Mutability ' + kmer)
            fig_dir_kmer = os.path.join(fig_dir, kmer)
            smart_makedirs(fig_dir_kmer)
            for chain in models:
                chain_grs = []
                fr_pdf = add_distr(betas_fr, chain, i, "FR")
                cdr_pdf = add_distr(betas_cdr, chain, i, "CDR")
                full_pdf = add_distr(betas_full, chain, i, "FULL")

                append_notnan(chain_grs, fr_pdf)
                append_notnan(chain_grs, cdr_pdf)
                append_notnan(chain_grs, full_pdf)
                append_notnan(grs_fr, fr_pdf)
                append_notnan(grs_cdr, cdr_pdf)
                append_notnan(grs_full, full_pdf)
                filename = os.path.join(fig_dir_kmer,
                                        "beta_" + chain.name + ".html")
                plot_fig(chain_grs, layout, filename)

            plot_fig(grs_fr, layout, os.path.join(fig_dir_kmer, 'beta_fr.html'))
            plot_fig(grs_cdr, layout, os.path.join(fig_dir_kmer, 'beta_cdr.html'))
            plot_fig(grs_full, layout, os.path.join(fig_dir_kmer, 'beta_full.html'))


    for strategy in MutationStrategies:
        str_figures_dir = os.path.join(figures_dir, strategy.name)
        smart_makedirs(str_figures_dir)
        print(strategy)
        plot_distr_str(models[strategy], str_figures_dir)


if __name__ == '__main__':
    input_config = read_config(parse_args().input)
    model_config = input_config.kmer_model_estimating
    models = read_models(model_config.outdir)
    figures_dir = os.path.join(model_config.figures_dir,
                               model_config.model_distribution_figs)
    plot_distr(models, figures_dir)
