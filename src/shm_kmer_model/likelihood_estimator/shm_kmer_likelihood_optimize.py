from shm_kmer_likelihood import ShmKmerLikelihood
from scipy.optimize import minimize
import numpy as np

class ShmKmerLikelihoodOptimizator:
    def __init__(self, shm_kmer_likelihood,
                 x0_beta=None, x0_dir=None,
                 bounds=None, method='L-BFGS-B'):
        if x0_beta is None or x0_dir is None:
            self.x0_beta, self.x0_dir = self.__get_start_point(shm_kmer_likelihood)
            # self.x0_beta, self.x0_dir = np.ones(2), np.ones(3)
        if bounds is None:
            self.bounds_beta=((0, None),) * 2
            self.bounds_dir=((0, None),) * shm_kmer_likelihood.mutated_sample.shape[1]
        self.method=method
        self.lkhd_beta = (lambda x: -shm_kmer_likelihood.likelihood_beta(x[0], x[1]))
        self.grad_beta = (lambda x: -shm_kmer_likelihood.gradient_beta(x[0], x[1]))
        self.lkhd_dir = (lambda x: -shm_kmer_likelihood.likelihood_dir(x))
        self.grad_dir = (lambda x: -shm_kmer_likelihood.gradient_dir(x))

    def __get_start_point(self, shm_kmer_likelihood, scale_beta=1, scale_dir=1):
        beta_shape1 = 1. - shm_kmer_likelihood.full_sample[:,shm_kmer_likelihood.nonmutated_ind]  / \
                          np.sum(shm_kmer_likelihood.full_sample, axis=1, dtype=float)
        beta_shape1 = np.mean(beta_shape1[~np.isnan(beta_shape1)])
        # beta_shape1 = np.sum(shm_kmer_likelihood.full_sample[:,shm_kmer_likelihood.nonmutated_ind])  / \
        #               np.sum(shm_kmer_likelihood.full_sample, dtype=float)
        beta_shape2 = 1. - beta_shape1
        beta_shape1 *= scale_beta
        beta_shape2 *= scale_beta

        scaled_mutated_sample = (shm_kmer_likelihood.mutated_sample.T / \
                                 np.sum(shm_kmer_likelihood.mutated_sample, axis=1, dtype=float)).T
        dir_lambda = np.mean(scaled_mutated_sample, axis=0, dtype=float)
        dir_lambda *= scale_dir
        return np.array([beta_shape1, beta_shape2]), dir_lambda
        # return np.array([beta_shape1, beta_shape2]), np.ones(3)

    def maximize(self):
        minimize_result_beta = minimize(fun=self.lkhd_beta,
                                        x0=self.x0_beta,
                                        bounds=self.bounds_beta,
                                        method=self.method,
                                        jac=self.grad_beta)
        minimize_result_dir = minimize(fun=self.lkhd_dir,
                                       x0=self.x0_dir,
                                       bounds=self.bounds_dir,
                                       method=self.method,
                                       jac=self.grad_dir)
        return self.x0_beta, self.x0_dir, minimize_result_beta, minimize_result_dir
