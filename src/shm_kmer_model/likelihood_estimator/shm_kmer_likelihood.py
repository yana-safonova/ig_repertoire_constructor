from scipy.special import gammaln, betaln, psi
from scipy.optimize import check_grad
from math_special_functions import multibetaln, dibetaln, dimultibetaln
import numpy as np

class ShmKmerLikelihood:
    def __init__(self, sample, nonmutated_ind, check_gradient=False, number_of_tests=None):
        self.full_sample = sample
        self.mutated_sample = np.delete(sample, nonmutated_ind, 1)
        self.nonmutated_ind = nonmutated_ind

        self.n_all = np.sum(self.full_sample, axis=1)
        self.n_mutated = np.sum(self.mutated_sample, axis=1)

        if check_gradient:
            if number_of_tests is None:
                print("Supply number of tests for check_gradient.")
            else:
                self.__check_gradient(number_of_tests)

    def likelihood(self, beta_shape1, beta_shape2, dir_lambda):
        result_p11 = betaln(beta_shape1 + self.n_mutated,
                            beta_shape2 + self.n_all - self.n_mutated)
        result_p12 = betaln(beta_shape1, beta_shape2)
        
        result_p21 = multibetaln(self.mutated_sample + dir_lambda, 1)
        result_p22 = multibetaln(dir_lambda)
            
        return np.sum(result_p11 - result_p12 + result_p21 - result_p22)

    def gradient(self, beta_shape1, beta_shape2, dir_lambda):
        N = self.full_sample.shape[0]
        result_p11 = np.sum(dibetaln(beta_shape1 + self.n_mutated,
                           beta_shape2 + self.n_all - self.n_mutated), 0)
        result_p12 = dibetaln(beta_shape1, beta_shape2).flatten() * N

        result_p21 = np.sum(dimultibetaln(self.mutated_sample + dir_lambda), 0)
        result_p22 = dimultibetaln(dir_lambda) * N

        return np.concatenate((result_p11 - result_p12,
                               result_p21 - result_p22))

    def __check_gradient(self, number_of_tests=None):
        cglikelihood = (lambda x: self.likelihood(beta_shape1=x[0],
                                                  beta_shape2=x[1],
                                                  dir_lambda=x[2:]))
        cggrad = (lambda x: self.gradient(beta_shape1=x[0],
                                          beta_shape2=x[1],
                                          dir_lambda=x[2:]))
        if number_of_tests is not None:
            for i in xrange(number_of_tests):
                point = np.random.exponential(size=5, scale=10)
                print(check_grad(cglikelihood, cggrad, x0=point))
        if x is not None:
            print(check_grad(cglikelihood, cggrad, x))
