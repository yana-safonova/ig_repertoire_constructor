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


from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from sample_reader.standard_samples import concatenate_kmer_matrices_all_data
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies
import numpy as np

matrices = concatenate_kmer_matrices_all_data()
matrices = matrices[MutationStrategies.Trivial][Chains.IGH]

filtered = filter_by_coverage(matrices, threshold_function=np.min)


from mutability_diversity.calculate_mutability_diversity import *
mut_fr = calculate_mutability_fr(matrices)
mut_cdr = calculate_mutability_cdr(matrices)
subst = calculate_substitution(matrices)
# print(np.sum(np.isnan(calculate_mutability_fr(matrices)), axis=1))

