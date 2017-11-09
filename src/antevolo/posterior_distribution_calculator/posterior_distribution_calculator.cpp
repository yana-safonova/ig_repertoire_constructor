//
// Created by Andrew Bzikadze on 3/10/17.
//

#include <unordered_set>

#include "posterior_distribution_calculator.hpp"
#include "kmer_utils/kmer_utils.hpp"

#include <iostream>

using namespace shm_kmer_matrix_estimator;

namespace antevolo {

ShmModel PosteriorDistributionCalculator::calculate(const ShmModel &prior,
                                                    const KmerMatrix &fr_matrix,
                                                    const KmerMatrix &cdr_matrix) const {
    ShmModel posterior(prior);
    ShmModel::ModelParametersVector posterior_beta_fr_params(prior.beta_fr_params());
    ShmModel::ModelParametersVector posterior_beta_cdr_params(prior.beta_cdr_params());
    ShmModel::ModelParametersVector posterior_beta_full_params(prior.beta_full_params());
    ShmModel::ModelParametersVector posterior_dirichlet_params(prior.dirichlet_params());

    ShmModel::ModelParametersVector posterior_start_point_beta_fr_params(prior.start_point_beta_fr_params());
    ShmModel::ModelParametersVector posterior_start_point_beta_cdr_params(prior.start_point_beta_cdr_params());
    ShmModel::ModelParametersVector posterior_start_point_beta_full_params(prior.start_point_beta_full_params());
    ShmModel::ModelParametersVector posterior_start_point_dirichlet_params(prior.start_point_dirichlet_params());

    std::vector<std::string> kmers(KmerUtils::GenerateAllKmersFixedLength(fr_matrix.kmer_len()));
    size_t central_index = fr_matrix.kmer_len() / 2;
    for (size_t i = 0; i < posterior.size(); ++i) {
        auto &kmer = kmers[i];

        auto &fr_sample = fr_matrix[i];
        auto &cdr_sample = cdr_matrix[i];
        KmerMatrixRowType full_sample;
        std::transform(fr_sample.cbegin(), fr_sample.cend(),
                       cdr_sample.cbegin(), full_sample.begin(),
                       std::plus<unsigned int>());

        auto &posterior_beta_fr_param = posterior_beta_fr_params[i];
        auto &posterior_beta_cdr_param = posterior_beta_cdr_params[i];
        auto &posterior_beta_full_param = posterior_beta_full_params[i];
        auto &posterior_dirichlet_param = posterior_dirichlet_params[i];

        auto &posterior_start_point_beta_fr_param = posterior_start_point_beta_fr_params[i];
        auto &posterior_start_point_beta_cdr_param = posterior_start_point_beta_cdr_params[i];
        auto &posterior_start_point_beta_full_param = posterior_start_point_beta_full_params[i];
        auto &posterior_start_point_dirichlet_param = posterior_start_point_dirichlet_params[i];

        size_t index_central_nucleotide = KmerUtils::GetIndexNucl(kmer[central_index]);
        auto update_beta_param = [&index_central_nucleotide]
            (decltype(fr_sample) sample,
             decltype(posterior_beta_fr_param)& beta_param) {
          size_t number_of_trials = sample[0] + sample[1] + sample[2] + sample[3];
          size_t number_of_successes = number_of_trials - sample[index_central_nucleotide];
          beta_param[0] += static_cast<double>(number_of_successes);
          beta_param[1] += static_cast<double>(number_of_trials - number_of_successes);
        };

        std::unordered_set<size_t> mutated_indices{0, 1, 2, 3};
        mutated_indices.erase(index_central_nucleotide);
        auto update_direchlet_param = [&mutated_indices, &full_sample]
            (decltype(posterior_dirichlet_param) dir_param) {
          size_t j = 0;
          for (const auto mutated_index : mutated_indices) {
              dir_param.at(j++) += full_sample[mutated_index];
          }
        };

        update_beta_param(fr_sample, posterior_beta_fr_param);
        update_beta_param(cdr_sample, posterior_beta_cdr_param);
        update_beta_param(full_sample, posterior_beta_full_param);
        update_direchlet_param(posterior_dirichlet_param);

        update_beta_param(fr_sample, posterior_start_point_beta_fr_param);
        update_beta_param(cdr_sample, posterior_start_point_beta_cdr_param);
        update_beta_param(full_sample, posterior_start_point_beta_full_param);
        update_direchlet_param(posterior_start_point_dirichlet_param);
    }

    posterior.set_beta_fr_params(posterior_beta_fr_params);
    posterior.set_beta_cdr_params(posterior_beta_cdr_params);
    posterior.set_beta_full_params(posterior_beta_full_params);
    posterior.set_dirichlet_params(posterior_dirichlet_params);

    posterior.set_start_point_beta_fr_params(posterior_start_point_beta_fr_params);
    posterior.set_start_point_beta_cdr_params(posterior_start_point_beta_cdr_params);
    posterior.set_start_point_beta_full_params(posterior_start_point_beta_full_params);
    posterior.set_start_point_dirichlet_params(posterior_start_point_dirichlet_params);

    return posterior;
}

} // End namespace antevolo
