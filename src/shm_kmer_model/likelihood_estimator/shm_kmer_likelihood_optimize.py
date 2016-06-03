from shm_kmer_likelihood import ShmKmerLikelihood
from scipy.optimize import minimize
import numpy as np

class ShmKmerLikelihoodOptimizator:
    def __init__(self, shm_kmer_likelihood, x0=None, bounds=None, method='L-BFGS-B'):
        if x0 is None:
            self.x0 = self.__get_start_point(shm_kmer_likelihood)
        if bounds is None:
            self.bounds=((0, None),) * (2 + shm_kmer_likelihood.mutated_sample.shape[1])
        self.method=method
        self.lkhd = (lambda x: -shm_kmer_likelihood.likelihood(beta_shape1=x[0],
                                                               beta_shape2=x[1],
                                                               dir_lambda=x[2:]))
        self.grad = (lambda x: -shm_kmer_likelihood.gradient(beta_shape1=x[0],
                                                             beta_shape2=x[1],
                                                             dir_lambda=x[2:]))
        
    def __get_start_point(self, shm_kmer_likelihood, scale_beta=1, scale_dir=1):
        beta_shape1 = np.mean(1 - shm_kmer_likelihood.full_sample[:,shm_kmer_likelihood.nonmutated_ind] / \
                              np.sum(shm_kmer_likelihood.full_sample, axis=1, dtype=float))
        beta_shape2 = 1 - beta_shape1
        beta_shape1 *= scale_beta
        beta_shape2 *= scale_beta
        
        scaled_mutated_sample = (shm_kmer_likelihood.mutated_sample.T / \
                                 np.sum(shm_kmer_likelihood.mutated_sample, axis=1, dtype=float)).T
        dir_lambda = np.mean(scaled_mutated_sample, axis=0, dtype=float)
        dir_lambda *= scale_dir
        return np.concatenate([np.array([beta_shape1, beta_shape2]), dir_lambda])
    
    def maximize(self):
        minimize_result = minimize(fun=self.lkhd,
                                   x0=self.x0,
                                   bounds=self.bounds,
                                   method=self.method,
                                   jac=self.grad)
        return minimize_result
