//
// Created by Andrew Bzikadze on 3/9/17.
//

#pragma once

#include <vector>
#include <fstream>

#include "kmer_utils/kmer_indexed_vector.hpp"

namespace antevolo {

class ShmModel {
private:
    using ModelParametersVector = shm_kmer_matrix_estimator::KmerIndexedVector<std::vector<double>>;

    ModelParametersVector beta_fr_params_;
    ModelParametersVector beta_cdr_params_;
    ModelParametersVector beta_full_params_;

    ModelParametersVector dirichlet_params_;

    ModelParametersVector start_point_beta_fr_params_;
    ModelParametersVector start_point_beta_cdr_params_;
    ModelParametersVector start_point_beta_full_params_;

    ModelParametersVector start_point_dirichlet_params_;

    using SuccessMLEOptimazation = shm_kmer_matrix_estimator::KmerIndexedVector<int>;
    SuccessMLEOptimazation beta_fr_success_mle_;
    SuccessMLEOptimazation beta_cdr_success_mle_;
    SuccessMLEOptimazation beta_full_success_mle_;
    SuccessMLEOptimazation dirichlet_success_mle_;

public:
    ShmModel() = delete;

    ShmModel(const ShmModel &) = default;
    ShmModel(ShmModel &&)      = default;
    ShmModel &operator=(const ShmModel &) = default;
    ShmModel &operator=(ShmModel &&)      = default;

    explicit ShmModel(const std::string& filename);

    virtual ~ShmModel() = default;

    const ModelParametersVector& beta_fr_params()   const { return beta_fr_params_; }
    const ModelParametersVector& beta_cdr_params()  const { return beta_cdr_params_; }
    const ModelParametersVector& beta_full_params() const { return beta_full_params_; }
    const ModelParametersVector& dirichlet_params() const { return dirichlet_params_; }

    const ModelParametersVector& start_point_beta_fr_params()   const { return start_point_beta_fr_params_; }
    const ModelParametersVector& start_point_beta_cdr_params()  const { return start_point_beta_cdr_params_; }
    const ModelParametersVector& start_point_beta_full_params() const { return start_point_beta_full_params_; }
    const ModelParametersVector& start_point_dirichlet_params() const { return start_point_dirichlet_params_; }

    const SuccessMLEOptimazation& beta_fr_success_mle()   const { return beta_fr_success_mle_; }
    const SuccessMLEOptimazation& beta_cdr_success_mle()  const { return beta_cdr_success_mle_; }
    const SuccessMLEOptimazation& beta_full_success_mle() const { return beta_full_success_mle_; }
    const SuccessMLEOptimazation& dirichlet_success_mle() const { return dirichlet_success_mle_; }

    size_t size() const;
};

} // End namespace antevolo
