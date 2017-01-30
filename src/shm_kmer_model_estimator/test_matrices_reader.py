from sample_reader.kmer_matrices_reader import KmerMatricesReader
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


reader = KmerMatricesReader()
matrices = reader.read("age")
print(matrices[MutationStrategies.Trivial][Chains.IGH].matrices[:, :, 0])

# from chains.chains import Chains
