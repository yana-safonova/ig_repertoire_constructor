//
// Created by Andrew Bzikadze on 5/15/16.
//

#pragma once
#include "shm_config.hpp"

namespace shm_kmer_model_estimator {
class SHMkmerModelEstimator {
    const shm_config::io_params &io_params_;
    const shm_config::alignment_checker_params &alignment_checker_params_;
    const shm_config::alignment_cropper_params &alignment_cropper_params_;
    const shm_config::mismatch_finder_params &mismatch_finder_params_;

public:
    SHMkmerModelEstimator(const shm_config::io_params &io_params,
                          const shm_config::alignment_checker_params &alignment_checker_params,
                          const shm_config::alignment_cropper_params &alignment_cropper_params,
                          const shm_config::mismatch_finder_params &mismatch_finder_params) :
        io_params_(io_params),
        alignment_checker_params_(alignment_checker_params),
        alignment_cropper_params_(alignment_cropper_params),
        mismatch_finder_params_(mismatch_finder_params) { }

    int Run();
};
}
