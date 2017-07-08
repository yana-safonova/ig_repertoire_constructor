//
// Created by Andrew Bzikadze on 5/18/16.
//

#pragma once

#include <vector>
#include <tuple>

#include "shm_kmer_matrix_estimator_config.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "alignment_checker/abstract_alignment_checker.hpp"
#include "alignment_checker/no_gaps_alignment_checker.hpp"
#include "alignment_cropper/abstract_alignment_cropper.hpp"
#include "alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp"

namespace shm_kmer_matrix_estimator {

class AlignmentReader {
private:
    std::string alignments_filename_;
    std::string cdr_details_filename_;
    AbstractAlignmentCheckerPtr alignment_checker_ptr_;
    AbstractAlignmentCropperPtr alignment_cropper_ptr_;

private:
    const size_t has_stop_codon_pos = 5;

public:
    AlignmentReader(const std::string &alignments_filename,
                    const std::string &cdr_details_filename,
                    const shm_kmer_matrix_estimator_config::alignment_checker_params &alignment_checker_params,
                    const shm_kmer_matrix_estimator_config::alignment_cropper_params &alignment_cropper_params);

    VectorEvolutionaryEdgeAlignments read_alignments() const;
};

} // End namespace shm_kmer_matrix_estimator
