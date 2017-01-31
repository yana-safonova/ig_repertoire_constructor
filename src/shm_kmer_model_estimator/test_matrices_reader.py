from sample_reader.kmer_matrices_reader import KmerMatricesReader
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies
from kmer_utilities.kmer_utilities import *


reader = KmerMatricesReader()
matrices = reader.read("age")
# print(matrices[MutationStrategies.Trivial][Chains.IGH].matrices[:, :, 0])
# print(matrices[MutationStrategies.Trivial][Chains.IGH][0:1].shape)

# from chains.chains import Chains

kmers = kmer_names()
for ind, kmer in enumerate(kmers):
    assert ind == kmer_index(kmer)

from sample_reader.standard_samples import concatenate_kmer_matrices_all_data

matrices = concatenate_kmer_matrices_all_data()
# print(matrices[MutationStrategies.Trivial][Chains.IGH]["AAAAA"])


#from config.config import config
#print(config.input_data)

from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood

test_sample = matrices[MutationStrategies.Trivial][Chains.IGH]["AAAAA"]
lklh = ShmKmerLikelihood(test_sample, check_gradient=True, number_of_tests=2)
