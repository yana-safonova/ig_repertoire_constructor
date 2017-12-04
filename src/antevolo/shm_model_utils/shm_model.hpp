//
// Created by Andrew Bzikadze on 3/9/17.
//

#pragma once

#include <vector>
#include <fstream>

#include "kmer_utils/kmer_indexed_vector.hpp"

#include "annotation_utils/annotated_clone.hpp"

namespace antevolo {

class ShmModel {
public:
    using ModelParametersVector = shm_kmer_matrix_estimator::KmerIndexedVector<std::vector<double>>;

private:
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
    const ModelParametersVector& beta_params(const annotation_utils::StructuralRegion&) const;
    const ModelParametersVector& dirichlet_params() const { return dirichlet_params_; }

    void set_beta_fr_params(const ModelParametersVector& beta_fr_params) { beta_fr_params_ = beta_fr_params; }
    void set_beta_cdr_params(const ModelParametersVector& beta_cdr_params) { beta_cdr_params_ = beta_cdr_params; }
    void set_beta_full_params(const ModelParametersVector& beta_full_params) { beta_full_params_ = beta_full_params; }
    void set_dirichlet_params(const ModelParametersVector& dirichlet_params) { dirichlet_params_ = dirichlet_params; }

    const ModelParametersVector& start_point_beta_fr_params()   const { return start_point_beta_fr_params_; }
    const ModelParametersVector& start_point_beta_cdr_params()  const { return start_point_beta_cdr_params_; }
    const ModelParametersVector& start_point_beta_full_params() const { return start_point_beta_full_params_; }
    const ModelParametersVector& start_point_beta_params(const annotation_utils::StructuralRegion& region) const;
    const ModelParametersVector& start_point_dirichlet_params() const { return start_point_dirichlet_params_; }

    void set_start_point_beta_fr_params(const ModelParametersVector& sp_fr) { start_point_beta_fr_params_ = sp_fr; }
    void set_start_point_beta_cdr_params(const ModelParametersVector& sp_cdr) { start_point_beta_cdr_params_ = sp_cdr; }
    void set_start_point_beta_full_params(const ModelParametersVector& sp_full) { start_point_beta_full_params_ = sp_full; }
    void set_start_point_dirichlet_params(const ModelParametersVector& sp_dir) { start_point_dirichlet_params_ = sp_dir; }

    const SuccessMLEOptimazation& beta_fr_success_mle()   const { return beta_fr_success_mle_; }
    const SuccessMLEOptimazation& beta_cdr_success_mle()  const { return beta_cdr_success_mle_; }
    const SuccessMLEOptimazation& beta_full_success_mle() const { return beta_full_success_mle_; }
    const SuccessMLEOptimazation& beta_success_mle(const annotation_utils::StructuralRegion&) const;
    const SuccessMLEOptimazation& dirichlet_success_mle() const { return dirichlet_success_mle_; }

    size_t size() const { return beta_fr_params_.size(); }
    size_t kmer_len() const { return beta_fr_params_.kmer_len(); }

    double beta_fr_expectation(const std::string& kmer) const;
    double beta_cdr_expectation(const std::string& kmer) const;
    double beta_full_expectation(const std::string& kmer) const;
    double beta_expectation(const std::string&,
                            const annotation_utils::StructuralRegion&) const;
    double dirichlet_expectation(const std::string& kmer, const char nucl) const;

    double likelihood_kmer_nucl(const std::string& kmer,
                                const char nucl,
                                const annotation_utils::StructuralRegion& region =
                                    annotation_utils::StructuralRegion::AnyRegion) const;
    double loglikelihood_kmer_nucl(const std::string& kmer,
                                   const char nucl,
                                   const annotation_utils::StructuralRegion& region =
                                       annotation_utils::StructuralRegion::AnyRegion) const;
};

std::ostream& operator<<(std::ostream& os, const ShmModel& obj);

} // End namespace antevolo
