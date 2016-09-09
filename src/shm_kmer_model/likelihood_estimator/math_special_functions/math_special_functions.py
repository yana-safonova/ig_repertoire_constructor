""" This module provides special routines for beta and dir lkhd calculation """
import numpy as np
from scipy.special import gammaln, psi


def multibetaln(x, axis=None):
    return np.sum(gammaln(x), axis=axis) - gammaln(np.sum(x, axis=axis))


def dibetaln(x, y):
    return (np.vstack((psi(x), psi(y))) - psi(x + y)).T


def dimultibetaln(x):
    if len(x.shape) == 2:
        second_part = psi(np.sum(x, axis=1))
    else:
        second_part = psi(np.sum(x))
    return (psi(x).T - second_part).T
