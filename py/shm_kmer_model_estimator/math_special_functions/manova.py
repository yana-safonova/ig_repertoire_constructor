import numpy as np
import scipy.stats as stats


class ManovaOutput(object):
    def __init__(self, pvalue, f_statistics):
        self.f_statistics = f_statistics
        self.pvalue = pvalue


def manova_oneway(*args):
    # args = [np.asarray(arg, dtype=float) for arg in args]
    lens = np.array([len(arg) for arg in args])

    all_data = np.concatenate(args)
    e_matrix = sum([np.cov(arg.T) * (len(arg) - 1) for arg in args])
    full_cov = np.cov(all_data.T) * (len(all_data) - 1)

    means = np.array([np.mean(arg, axis=0) for arg in args])
    big_mean = np.mean(all_data, axis=0)

    h_matrix = np.zeros((len(big_mean), len(big_mean)))
    for mean, single_len in zip(means, lens):
        centr_mean = mean - big_mean
        h_matrix += np.outer(centr_mean, centr_mean) * single_len

    return e_matrix, full_cov


def manova_wilks(*args):
    if len(args) == 0:
        return
    e, full_cov = manova_oneway(*args)
    lw = np.linalg.det(e) / np.linalg.det(full_cov)
    n = sum([len(arg) for arg in args])
    k = len(args)
    p = args[0].shape[1]

    nu_h = k - 1
    nu_e = n - k
    f = nu_e - 0.5 * (p - nu_h + 1)
    g = -1 + 0.5 * p * nu_h
    if p**2 + nu_h**2 - 5 > 0:
        t = np.sqrt((p**2 * nu_h**2 - 4.) / (p**2 + nu_h**2 - 5.))
    else:
        t = 1

    f_stat = (f * t - g) * (1 - lw**(1/t)) / (p * nu_h * lw**(1/t))
    pvalue = 1 - stats.f.cdf(f_stat, p * nu_h, f * t - g)
    return ManovaOutput(pvalue, f_stat)


# a = np.array([[7.0, 17.0, 19.7],
#               [7.3, 17.2, 20.3],
#               [8.0, 19.3, 22.6],
#               [8.1, 19.8, 23.7],
#               [7.9, 18.4, 22.0]])
# b = np.array([[7.3, 17.4, 22.5],
#               [7.7, 19.8, 24.9],
#               [8.2, 20.2, 26.1],
#               [8.3, 22.6, 27.5],
#               [6.4, 23.4, 28.1]])
# print(manova_wilks((a, b)))
