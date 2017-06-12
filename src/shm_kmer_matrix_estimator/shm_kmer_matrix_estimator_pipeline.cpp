//
// Created by Andrew Bzikadze on 5/15/16.
//

#include <algorithm>

#include "logger/logger.hpp"
#include "shm_kmer_matrix_estimator_pipeline.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "germline_alignment_reader/alignment_reader.hpp"
#include "shm_kmer_matrix_estimator/shm_kmer_matrix_estimator.hpp"
#include "kmer_matrix_exporter/kmer_matrix_exporter.hpp"

namespace shm_kmer_matrix_estimator {

int SHMkmerMatrixEstimatorPipeline::Run() const {
    const std::string boarder("=============");
    INFO(boarder << " SHM k-mer Model Calculator starts " << boarder);

    INFO(boarder << " Reading alignments starts " << boarder);
    INFO(std::string("input v_alignments filename: ") << io_params_.input.v_alignments);
    INFO(std::string("input cdr_details filename: ") << io_params_.input.cdr_details);

    INFO(std::string("Strategy for checking: ") <<
        alignment_checker_params_.alignment_checker_method_names[
            static_cast<size_t> (alignment_checker_params_.alignment_checker_method)]);
    INFO(std::string("Strategy for cropping: ") <<
        alignment_cropper_params_.alignment_cropper_method_names[
            static_cast<size_t> (alignment_cropper_params_.alignment_cropper_method)]);

    AlignmentReader germline_alignment_reader(io_params_.input.v_alignments,
                                              io_params_.input.cdr_details,
                                              alignment_checker_params_,
                                              alignment_cropper_params_);
    VectorEvolutionaryEdgeAlignments alignments(germline_alignment_reader.read_alignments());
    INFO(boarder << " Reading alignments finishes " << boarder);
    INFO(boarder << " Read " << alignments.size() << " sequences " << boarder);

    INFO(boarder << " Estimating statistics starts " << boarder);
    INFO(std::string("Strategy for mutations: ") <<
        mutations_strategy_params_.mutation_strategy_method_names[
            static_cast<size_t>(mutations_strategy_params_.mutations_strategy_method)]);

    ShmKmerMatrixEstimator statistics_estimator(mutations_strategy_params_,
                                             alignment_checker_params_,
                                             alignment_cropper_params_);
    KmerMatrix statistics_fr, statistics_cdr;
    std::tie(statistics_fr, statistics_cdr) = statistics_estimator.calculate_mutation_statistics(alignments);
    INFO(boarder << " Estimating statistics finishes" << boarder);

    INFO(boarder << " Exporting statistics starts " << boarder);
    INFO(std::string("Output filename for FR: ") << io_params_.output.output_filename_fr);
    INFO(std::string("Output filename for CDR: ") << io_params_.output.output_filename_cdr);
    KmerMatrixExporter statistics_exporter;
    statistics_exporter.export_statistics(io_params_.output.output_filename_fr, statistics_fr);
    statistics_exporter.export_statistics(io_params_.output.output_filename_cdr, statistics_cdr);
    INFO(boarder << " Exporting statistics finishes " << boarder);

    return 0;
}

} // End namespace shm_kmer_matrix_estimator
