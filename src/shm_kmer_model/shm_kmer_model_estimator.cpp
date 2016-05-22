//
// Created by Andrew Bzikadze on 5/15/16.
//

#include "logger/logger.hpp"
#include "shm_kmer_model_estimator.hpp"

#include "gene_alignment.hpp"
#include "alignment_reader.hpp"

int shm_kmer_model_estimator::SHMkmerModelEstimator::Run() {
    ns_alignment_reader::AlignmentReader alignment_reader(io_params_.input.input_filename,
                                                          alignment_checker_params_,
                                                          alignment_cropper_params_);
    ns_gene_alignment::VectorReadGermlinePairs alignments(alignment_reader.read_alignments());
    // for (auto& alignment : alignments) {
    //     INFO(alignment.first);
    //     INFO(alignment.second);
    // }
    return 0;
}