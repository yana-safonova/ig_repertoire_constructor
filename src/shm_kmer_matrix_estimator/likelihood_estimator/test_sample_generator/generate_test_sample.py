import numpy as np


def generate_sample(sample_size, number_samples,
                    real_beta_shape1, real_beta_shape2, real_dir_lambda,
                    nonmutated_ind):
    mutation_bin_prob = np.random.beta(real_beta_shape1, real_beta_shape2,
                                       size=number_samples)
    is_mutated = np.random.binomial(sample_size, p=mutation_bin_prob,
                                    size=number_samples)

    mutation_dir_probs = np.random.dirichlet(real_dir_lambda,
                                             size=number_samples)
    sample = []
    for n, pvals in zip(is_mutated, mutation_dir_probs):
        sample.append(np.random.multinomial(n=n, pvals=pvals, size=1))
    sample = np.array(sample).reshape((number_samples, len(real_dir_lambda)))
    final_sample = np.insert(sample, nonmutated_ind,
                             sample_size - is_mutated, axis=1)
    return final_sample
