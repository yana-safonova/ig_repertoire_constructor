from __future__ import division
from scipy.optimize import minimize
import numpy as np


class ShmKmerLikelihoodOptimizator:
    """ This class is for lkhd optimization.
    Input: lklh itself. """
    def __init__(self, lklh,
                 x0_beta_fr=None, x0_beta_cdr=None,
                 x0_beta_full=None, x0_dir=None,
                 bounds=None, method='L-BFGS-B'):
        self.x0_beta_fr, self.x0_beta_cdr, self.x0_beta_full, self.x0_dir = \
            self.__get_start_point(lklh)
        if bounds is None:
            positive = ((0, None),)
            self.bounds_beta = positive * 2
            self.bounds_dir = positive * 3
        self.method = method

        self.lklh = lklh
        self.lkhd_beta_fr = (lambda x: -lklh.likelihood_beta_fr(*x))
        self.grad_beta_fr = (lambda x: -lklh.gradient_beta_fr(*x))
        self.lkhd_beta_cdr = (lambda x: -lklh.likelihood_beta_cdr(*x))
        self.grad_beta_cdr = (lambda x: -lklh.gradient_beta_cdr(*x))
        self.lkhd_beta_full = (lambda x: -lklh.likelihood_beta_full(*x))
        self.grad_beta_full = (lambda x: -lklh.gradient_beta_full(*x))
        self.lkhd_dir = (lambda x: -lklh.likelihood_dir(x))
        self.grad_dir = (lambda x: -lklh.gradient_dir(x))

    def __get_start_point(self, lklh, scale_beta=None, scale_dir=None):
        def get_beta_start_point(sample, mut, nmut, scale_beta=scale_beta):
            means = sample[:, mut] / (sample[:, mut] + sample[:, nmut])
            means = means[~np.isnan(means)]
            if len(means) == 0:
                return np.repeat(np.nan, 2)

            first_moment = np.nanmean(means)
            beta_shape = np.nanmean(means)

            if scale_beta is None:
                second_moment = np.nanmean(means**2)
                if np.allclose(second_moment, first_moment):
                    return np.repeat(np.nan, 2)
                scale_beta = \
                    (first_moment - second_moment) / \
                    (second_moment - first_moment**2)
                assert scale_beta > 0

            beta_shape = \
                scale_beta * np.array([first_moment, 1 - first_moment])
            return beta_shape

        beta_shape_fr = get_beta_start_point(sample=lklh.sample,
                                             mut=lklh.fr_mut_ind,
                                             nmut=lklh.fr_nonmut_ind)
        beta_shape_cdr = get_beta_start_point(sample=lklh.sample,
                                              mut=lklh.cdr_mut_ind,
                                              nmut=lklh.cdr_nonmut_ind)
        sample = np.copy(lklh.sample)
        sample[:, lklh.fr_mut_ind] = \
            np.sum(sample[:, [lklh.fr_mut_ind, lklh.cdr_mut_ind]], axis=1)
        sample[:, lklh.fr_nonmut_ind] = \
            np.sum(sample[:, [lklh.fr_nonmut_ind, lklh.cdr_nonmut_ind]], axis=1)
        beta_shape_full = get_beta_start_point(sample=sample,
                                               mut=lklh.fr_mut_ind,
                                               nmut=lklh.fr_nonmut_ind)

        scaled_mutated_sample = (lklh.mutated_sample.T /
                                 np.sum(lklh.mutated_sample, axis=1)).T
        dir_mean = np.nanmean(scaled_mutated_sample, axis=0)
        if scale_dir is None:
            if np.any(dir_mean == 1) or \
               np.any(dir_mean == 0) or \
               np.any(np.isnan(dir_mean)):
                scale_dir = np.nan
            else:
                dir_var = np.nanvar(scaled_mutated_sample, axis=0)
                scale_dir = dir_var / (dir_mean * (1 - dir_mean))
                scale_dir = np.nanmean(scale_dir)
                assert scale_dir > 0

        dir_lambda = scale_dir * dir_mean
        return beta_shape_fr, beta_shape_cdr, beta_shape_full, dir_lambda

    def maximize(self):
        def smart_minimize(*args, **kwargs):
            return minimize(*args, **kwargs)
        minimize_result_beta_fr = smart_minimize(fun=self.lkhd_beta_fr,
                                                 x0=self.x0_beta_fr,
                                                 bounds=self.bounds_beta,
                                                 method=self.method,
                                                 jac=self.grad_beta_fr)
        minimize_result_beta_cdr = smart_minimize(fun=self.lkhd_beta_cdr,
                                                  x0=self.x0_beta_cdr,
                                                  bounds=self.bounds_beta,
                                                  method=self.method,
                                                  jac=self.grad_beta_cdr)
        minimize_result_beta_full = smart_minimize(fun=self.lkhd_beta_full,
                                                  x0=self.x0_beta_full,
                                                  bounds=self.bounds_beta,
                                                  method=self.method,
                                                  jac=self.grad_beta_full)
        minimize_result_dir = smart_minimize(fun=self.lkhd_dir,
                                             x0=self.x0_dir,
                                             bounds=self.bounds_dir,
                                             method=self.method,
                                             jac=self.grad_dir)
        # if not minimize_result_dir.success:
        #     print self.lklh.sample[:, 4:]
        #     print self.x0_dir
        #     print self.lkhd_dir(self.x0_dir)
        #     print minimize_result_dir
        #     assert False
        return {'start_beta_fr': self.x0_beta_fr,
                'start_beta_cdr': self.x0_beta_cdr,
                'start_beta_full': self.x0_beta_full,
                'start_dir': self.x0_dir,
                'optim_res_beta_fr': minimize_result_beta_fr,
                'optim_res_beta_cdr': minimize_result_beta_cdr,
                'optim_res_beta_full': minimize_result_beta_full,
                'optim_res_dir': minimize_result_dir}
