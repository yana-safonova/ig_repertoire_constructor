//
// Created by Andrew Bzikadze on 5/18/16.
//

#pragma once

#include <vector>
#include <tuple>

#include "shm_config.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "alignment_checker/abstract_alignment_checker.hpp"
#include "alignment_checker/no_gaps_alignment_checker.hpp"
#include "alignment_cropper/abstract_alignment_cropper.hpp"
#include "alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp"

namespace ns_alignment_reader {

class AlignmentReader {
private:
    std::string alignments_filename_;
    std::string cdr_details_filename_;
    ns_abstract_alignment_checker::AbstractAlignmentCheckerPtr alignment_checker_ptr_;
    ns_abstract_alignment_cropper::AbstractAlignmentCropperPtr alignment_cropper_ptr_;

private:
    const size_t csv_cdr1_const = 9;

public:
    AlignmentReader(const std::string &alignments_filename,
                    const std::string &cdr_details_filename,
                    const shm_config::alignment_checker_params &alignment_checker_params,
                    const shm_config::alignment_cropper_params &alignment_cropper_params);

    ns_gene_alignment::VectorEvolutionaryEdgeAlignments read_alignments() const;
};
}
