# from sample_reader.kmer_matrices_reader import KmerMatricesReader
# from chains.chains import Chains
# from mutation_strategies.mutation_strategies import MutationStrategies
# from kmer_utilities.kmer_utilities import *
# 
# 
# reader = KmerMatricesReader()
# matrices = reader.read("age")
# # print(matrices[MutationStrategies.Trivial][Chains.IGH].matrices[:, :, 0])
# # print(matrices[MutationStrategies.Trivial][Chains.IGH][0:1].shape)
# 
# # from chains.chains import Chains
# 
# kmers = kmer_names()
# for ind, kmer in enumerate(kmers):
#     assert ind == kmer_index(kmer)
# 
# from sample_reader.standard_samples import concatenate_kmer_matrices_all_data
# 
# matrices = concatenate_kmer_matrices_all_data()
# # print(matrices[MutationStrategies.Trivial][Chains.IGH]["AAAAA"])
# # 
# # 
# # #from config.config import config
# # #print(config.input_data)
# # 
# # from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
# # 
# # test_sample = matrices[MutationStrategies.Trivial][Chains.IGH]["AAAAA"]
# # lklh = ShmKmerLikelihood(test_sample, check_gradient=True, number_of_tests=2)
# # 
# # from test_sample_generator.generate_test_sample import generate_sample
# # from shm_kmer_likelihood_optimize.shm_kmer_likelihood_optimize import ShmKmerLikelihoodOptimizator
# # from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
# # 
# # sample = generate_sample(100, 40, [1, 5], [1, 1], [20, 1, 10])
# # lklh = ShmKmerLikelihood(sample)
# # optimizator = ShmKmerLikelihoodOptimizator(lklh)
# # optim_res = optimizator.maximize()
# # print(optim_res['optim_result_beta_fr'].x)
# # print(optim_res['optim_result_beta_cdr'].x)
# # print(optim_res['optim_res_dir'].x)
# 
# 
# # from shm_kmer_model.yale_model import YaleSHM_Model, YaleMethod
# # 
# # yale_model = YaleSHM_Model()
# # print(yale_model.all_expect_mut_probs(YaleMethod.Inferred))
# 
# from shm_kmer_model_estimator.shm_kmer_model_estimator import ShmKmerModelEstimator
# 
# estimator = ShmKmerModelEstimator()
# 
# from chains.chains import Chains
# from mutation_strategies.mutation_strategies import MutationStrategies
# import pickle
# from config.config import config
# 
# print(config.output_csv_header)
# 
# models = estimator.estimate_models(matrices)
# for strategy in models:
#     for chain in models[strategy]:
#         model = models[strategy][chain]
#         model.dataset.to_csv('models/' + strategy.name + '_' + chain.name + '.csv', na_rep='nan')


# from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
# from sample_reader.standard_samples import concatenate_kmer_matrices_all_data
# from chains.chains import Chains
# from mutation_strategies.mutation_strategies import MutationStrategies
# import numpy as np
# 
# matrices = concatenate_kmer_matrices_all_data()
# matrices = matrices[MutationStrategies.Trivial][Chains.IGH]
# 
# filtered = filter_by_coverage(matrices, threshold_function=np.min)
# 
# 
# from mutability_diversity.calculate_mutability_diversity import *
# mut_fr = calculate_mutability_fr(matrices)
# mut_cdr = calculate_mutability_cdr(matrices)
# subst = calculate_substitution(matrices)
# print(np.sum(np.isnan(calculate_mutability_fr(matrices)), axis=1))


from spots.spots import coldspots, hotspots

from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies
from model_estimator import read_input_config
from mutability_diversity.calculate_mutability_diversity import calculate_mutability_fr, calculate_mutability_cdr, calculate_mutability_full
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


prefix_dir, input_data, model_est_params = read_input_config()
matrices = concatenate_kmer_freq_matrices(
    input_data=input_data,
    prefix_dir=prefix_dir,
    dir_data=model_est_params.kmer_matrices_dir,
    filename_fr=model_est_params.filename_fr,
    filename_cdr=model_est_params.filename_cdr)

def plot_mutability_boxplots(matrices, output_dir):
    hspots = set(hotspots())
    cspots = set(coldspots())

    def plot_by_strategy_chain(matrices, chain, strategy):
        matrices = filter_by_coverage(
            matrices, 
            threshold_function = lambda x, axis: np.min(x[:, [1, 3]], axis=axis),
            coverage_threshold=1000)
        mut_fr = calculate_mutability_fr(matrices)
        mut_cdr = calculate_mutability_cdr(matrices)
        mut_full = calculate_mutability_full(matrices)

        def boxplots(mut, kmer_names, fig_name):
            med = np.median(mut, axis=1)
            order = np.argsort(med)[::-1]
            kmer_names = np.array(kmer_names)
            kmer_names = kmer_names[order]
            mut = mut[order, :]
            mut = pd.DataFrame(mut.T, columns=kmer_names)

            colors = np.array(['black'] * len(kmer_names))
            for i, kmer in enumerate(kmer_names):
                if kmer in hspots:
                    colors[i] = 'red'
                elif kmer in cspots:
                    colors[i] = 'blue'
            g = sns.boxplot(mut)
            g.set_xticklabels(g.get_xticklabels(), rotation=90)
            [t.set_color(i) for (i,t) in zip(colors, g.xaxis.get_ticklabels())]
            [t.set_facecolor(i) for (i,t) in zip(colors, g.artists)]
            fig = g.get_figure()
            fig.savefig(fig_name, format='pdf')
            plt.close()

        kmer_names = matrices.kmer_names
        fig_name_fr = "output_%s_%s_fr.pdf" % (strategy.name, chain.name)
        fig_name_cdr = "output_%s_%s_cdr.pdf" % (strategy.name, chain.name)
        fig_name_full = "output_%s_%s_full.pdf" % (strategy.name, chain.name)
        boxplots(mut_fr, kmer_names, fig_name_fr)
        boxplots(mut_cdr, kmer_names, fig_name_cdr)
        boxplots(mut_full, kmer_names, fig_name_full)

    for strategy in MutationStrategies:
        for chain in Chains:
            print(strategy, chain)
            plot_by_strategy_chain(matrices[strategy][chain],
                                   chain, strategy)

plot_mutability_boxplots(matrices, "")
