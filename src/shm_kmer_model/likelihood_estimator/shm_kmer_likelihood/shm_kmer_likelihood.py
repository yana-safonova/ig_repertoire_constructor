from scipy.special import gammaln, betaln, psi
from scipy.optimize import check_grad
from math_special_functions.math_special_functions import multibetaln, dibetaln, dimultibetaln
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
                self.check_gradient(number_of_tests)

    def likelihood_beta(self, beta_shape1, beta_shape2):
        result_p11 = betaln(beta_shape1 + self.n_mutated,
                            beta_shape2 + self.n_all - self.n_mutated)
        result_p12 = betaln(beta_shape1, beta_shape2)

        return np.sum(result_p11 - result_p12)

    def likelihood_dir(self, dir_lambda):
        result_p21 = multibetaln(self.mutated_sample + dir_lambda, 1)
        result_p22 = multibetaln(dir_lambda)

        return np.sum(result_p21 - result_p22)


    def gradient_beta(self, beta_shape1, beta_shape2):
        N = self.full_sample.shape[0]
        result_p11 = np.sum(dibetaln(beta_shape1 + self.n_mutated,
                           beta_shape2 + self.n_all - self.n_mutated), 0)
        result_p12 = dibetaln(beta_shape1, beta_shape2).flatten() * N

        return result_p11 - result_p12

    def gradient_dir(self, dir_lambda):
        N = self.full_sample.shape[0]

        result_p21 = np.sum(dimultibetaln(self.mutated_sample + dir_lambda), 0)
        result_p22 = dimultibetaln(dir_lambda) * N

        return result_p21 - result_p22

    def check_gradient(self, number_of_tests=None):
        cglikelihood_p0 = (lambda x: self.likelihood_beta(x[0], x[1]))
        cggrad_p0 = (lambda x: self.gradient_beta(x[0], x[1]))

        cglikelihood_p1 = (lambda x: self.likelihood_dir(x))
        cggrad_p1 = (lambda x: self.gradient_dir(x))

        if number_of_tests is not None:
            for i in xrange(number_of_tests):
                point = np.random.exponential(size=5, scale=10)
                print(check_grad(cglikelihood_p0, cggrad_p0, x0=point[:2]))
                print(check_grad(cglikelihood_p1, cggrad_p1, x0=point[2:]))
