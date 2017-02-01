import numpy as np


def generate_sample(sample_size, n_samples,
                    beta_fr, beta_cdr, dir_lambda):
    def get_binomial(beta):
        bin_p = np.random.beta(*beta, size=n_samples)
        mut = np.random.binomial(sample_size, p=bin_p, size=n_samples)
        mut = mut[:, np.newaxis]
        nmut = sample_size - mut
        return np.hstack((mut, nmut))

    fr = get_binomial(beta_fr)
    cdr = get_binomial(beta_cdr)

    dir_p = np.random.dirichlet(dir_lambda, size=n_samples)
    sample = []
    for i in xrange(n_samples):
        nmut = fr[i, 1] + cdr[i, 1]
        s = np.random.multinomial(n=nmut, pvals=dir_p[i], size=1)
        sample.append(s)
    mut_sample = np.array(sample).reshape((n_samples, len(dir_lambda)))
    final_sample = np.concatenate((fr, cdr, mut_sample), axis=1)
    return final_sample
