from scipy.special import betaln
from scipy.optimize import check_grad
from math_special_functions.math_special_functions \
        import multibetaln, dibetaln, dimultibetaln
import numpy as np


class ShmKmerLikelihood:
    """ This class calculates likelihood for a kmer.
    sample is np.array of shape #n_datasets * 5
    No optimization is provided here. """
    def __init__(self, sample,
                check_gradient=False, number_of_tests=None):
        self.sample = sample
        self.fr_mut_ind, self.fr_nonmut_ind = 0, 1
        self.cdr_mut_ind, self.cdr_nonmut_ind = 2, 3
        self.mutated_sample = sample[:, (self.cdr_nonmut_ind + 1):]

        if check_gradient:
            if number_of_tests is None:
                raise ValueError("Supply number of tests for check_gradient.")
            else:
                self.check_gradient(number_of_tests)

    def __likelihood_beta(self, beta_shape1, beta_shape2, mut, nmut):
        result_p11 = betaln(beta_shape1 + mut,
                            beta_shape2 + nmut)
        result_p12 = betaln(beta_shape1, beta_shape2)

        return np.sum(result_p11 - result_p12)

    def likelihood_beta_fr(self, beta_shape1, beta_shape2):
        return self.__likelihood_beta(beta_shape1=beta_shape1,
                                      beta_shape2=beta_shape2,
                                      mut=self.sample[:, self.fr_mut_ind],
                                      nmut=self.sample[:, self.fr_nonmut_ind])

    def likelihood_beta_cdr(self, beta_shape1, beta_shape2):
        return self.__likelihood_beta(beta_shape1=beta_shape1,
                                      beta_shape2=beta_shape2,
                                      mut=self.sample[:, self.cdr_mut_ind],
                                      nmut=self.sample[:, self.cdr_nonmut_ind])

    def likelihood_beta_full(self, beta_shape1, beta_shape2):
        mut = np.sum(self.sample[:, [self.fr_mut_ind,
                                     self.cdr_mut_ind]], axis=1)
        nmut = np.sum(self.sample[:, [self.fr_nonmut_ind,
                                      self.cdr_nonmut_ind]], axis=1)
        return self.__likelihood_beta(beta_shape1=beta_shape1,
                                      beta_shape2=beta_shape2,
                                      mut=mut,
                                      nmut=nmut)

    def likelihood_dir(self, dir_lambda):
        result_p21 = multibetaln(self.mutated_sample + dir_lambda, 1)
        result_p22 = multibetaln(dir_lambda)

        return np.sum(result_p21 - result_p22)

    def __gradient_beta(self, beta_shape1, beta_shape2, mut, nmut):
        N = self.sample.shape[0]
        result_p11 = np.sum(dibetaln(beta_shape1 + mut,
                                     beta_shape2 + nmut), 0)
        result_p12 = dibetaln(beta_shape1, beta_shape2).flatten() * N

        return result_p11 - result_p12

    def gradient_beta_fr(self, beta_shape1, beta_shape2):
        return self.__gradient_beta(beta_shape1, beta_shape2,
                                    mut=self.sample[:, self.fr_mut_ind],
                                    nmut=self.sample[:, self.fr_nonmut_ind])

    def gradient_beta_cdr(self, beta_shape1, beta_shape2):
        return self.__gradient_beta(beta_shape1, beta_shape2,
                                    mut=self.sample[:, self.cdr_mut_ind],
                                    nmut=self.sample[:, self.cdr_nonmut_ind])

    def gradient_beta_full(self, beta_shape1, beta_shape2):
        mut = np.sum(self.sample[:, [self.fr_mut_ind,
                                     self.cdr_mut_ind]], axis=1)
        nmut = np.sum(self.sample[:, [self.fr_nonmut_ind,
                                      self.cdr_nonmut_ind]], axis=1)
        return self.__gradient_beta(beta_shape1=beta_shape1,
                                    beta_shape2=beta_shape2,
                                    mut=mut,
                                    nmut=nmut)

    def gradient_dir(self, dir_lambda):
        N = self.sample.shape[0]

        result_p21 = np.sum(dimultibetaln(self.mutated_sample + dir_lambda), 0)
        result_p22 = dimultibetaln(dir_lambda) * N

        return result_p21 - result_p22

    def check_gradient(self, number_of_tests, tol=0.1):
        cglikelihood_p0_fr = (lambda x: self.likelihood_beta_fr(*x))
        cggrad_p0_fr = (lambda x: self.gradient_beta_fr(*x))

        cglikelihood_p0_cdr = (lambda x: self.likelihood_beta_cdr(*x))
        cggrad_p0_cdr = (lambda x: self.gradient_beta_cdr(*x))

        cglikelihood_p0_full = (lambda x: self.likelihood_beta_full(*x))
        cggrad_p0_full = (lambda x: self.gradient_beta_full(*x))

        cglikelihood_p1 = (lambda x: self.likelihood_dir(x))
        cggrad_p1 = (lambda x: self.gradient_dir(x))

        for i in xrange(number_of_tests):
            point = np.random.exponential(size=7, scale=10)
            assert check_grad(cglikelihood_p0_full, cggrad_p0_full, x0=point[:2]) < tol
            assert check_grad(cglikelihood_p0_fr, cggrad_p0_fr, x0=point[:2]) < tol
            assert check_grad(cglikelihood_p0_cdr, cggrad_p0_cdr, x0=point[2:4]) < tol
            assert check_grad(cglikelihood_p1, cggrad_p1, x0=point[4:]) < tol
