from __future__ import division
from scipy.optimize import minimize
import numpy as np


class ShmKmerLikelihoodOptimizator:
    """ This class is for lkhd optimization.
    Input: lklh itself. """
    def __init__(self, lklh,
                 x0_beta_fr=None, x0_beta_cdr=None, x0_dir=None,
                 bounds=None, method='L-BFGS-B'):
        if x0_beta_fr is None or x0_beta_cdr is None or x0_dir is None:
            self.x0_beta_fr, self.x0_beta_cdr, self.x0_dir = \
                self.__get_start_point(lklh)
        # print(self.x0_beta_fr, self.x0_beta_cdr, self.x0_dir)
        if bounds is None:
            positive = ((0, None),)
            self.bounds_beta = positive * 2
            self.bounds_dir = positive * 3
        self.method = method

        self.lkhd_beta_fr = (lambda x: -lklh.likelihood_beta_fr(*x))
        self.grad_beta_fr = (lambda x: -lklh.gradient_beta_fr(*x))
        self.lkhd_beta_cdr = (lambda x: -lklh.likelihood_beta_cdr(*x))
        self.grad_beta_cdr = (lambda x: -lklh.gradient_beta_cdr(*x))
        self.lkhd_dir = (lambda x: -lklh.likelihood_dir(x))
        self.grad_dir = (lambda x: -lklh.gradient_dir(x))

    def __get_start_point(self, lklh, scale_beta=1, scale_dir=1):
        def get_beta_start_point(sample, mut, nmut, scale_beta=scale_beta):
            beta_shape = sample[:, mut] / (sample[:, mut] + sample[:, nmut])
            beta_shape = np.mean(beta_shape[~np.isnan(beta_shape)])
            beta_shape = scale_beta * np.array([beta_shape, 1 - beta_shape])
            return beta_shape

        beta_shape_fr = get_beta_start_point(sample=lklh.sample,
                                             mut=lklh.fr_mut_ind,
                                             nmut=lklh.fr_nonmut_ind)
        beta_shape_cdr = get_beta_start_point(sample=lklh.sample,
                                              mut=lklh.cdr_mut_ind,
                                              nmut=lklh.cdr_nonmut_ind)

        scaled_mutated_sample = (lklh.mutated_sample.T /
                                 np.sum(lklh.mutated_sample, axis=1)).T
        dir_lambda = np.mean(scaled_mutated_sample, axis=0)
        dir_lambda *= scale_dir
        return beta_shape_fr, beta_shape_cdr, dir_lambda

    def maximize(self):
        minimize_result_beta_fr = minimize(fun=self.lkhd_beta_fr,
                                           x0=self.x0_beta_fr,
                                           bounds=self.bounds_beta,
                                           method=self.method,
                                           jac=self.grad_beta_fr)
        minimize_result_beta_cdr = minimize(fun=self.lkhd_beta_cdr,
                                           x0=self.x0_beta_cdr,
                                           bounds=self.bounds_beta,
                                           method=self.method,
                                           jac=self.grad_beta_cdr)
        minimize_result_dir = minimize(fun=self.lkhd_dir,
                                       x0=self.x0_dir,
                                       bounds=self.bounds_dir,
                                       method=self.method,
                                       jac=self.grad_dir)
        return {'start_beta_fr': self.x0_beta_fr,
                'start_beta_cdr': self.x0_beta_cdr,
                'start_dir': self.x0_dir,
                'optim_res_beta_fr': minimize_result_beta_fr,
                'optim_res_beta_cdr': minimize_result_beta_cdr,
                'optim_res_dir': minimize_result_dir}
