class LikelihoodOptimizator:
    def __init__(self, sample, nonmutated_ind):
        self.sample = sample
        self.nonmutated_ind = nonmutated_ind

    def likelihood(beta_shape1, beta_shape2, dir_lambda):
        from scipy.special import betaln
        sum_sample_lambda = self.sample + dir_lambda
        result_p1 = np.sum(scipy.special(sum_sample_lambda), axis=1)
        result_p1 -= betaln(np.sum(sum_sample_lambda, axis=1))
        result_p1 = np.sum(result, axis=0)

        n_all = np.sum(self.sample, axis=0)
        n_mutated = np.sum(np.delete(self.sample, nonmutated_ind, 0), axis=0)
        result_p2 = betaln(beta_shape1 + n_mutated)
        result_p2 += betaln(beta_shape2 + n_all - n_mutated)
        result_p2 -= betaln(beta_shape1 + beta_shape2 + n_all)
        result_p2 = np.sum(result_p2)

        result_p3 = np.sum(betaln(dir_lambda))
        result_p3 -= betaln(np.sum(dir_lambda))
        result_p3 += betaln(beta_shape1) + betaln(beta_shape2)
        result_p3 -= betaln(beta_shape1 + beta_shape2)
        result_p3 *= -self.sample.shape[1]

        return result_p1 + result_p2 + result_p3
