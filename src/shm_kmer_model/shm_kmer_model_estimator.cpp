//
// Created by Andrew Bzikadze on 5/15/16.
//

#include <algorithm>

#include "logger/logger.hpp"
#include "shm_kmer_model_estimator.hpp"
#include "gene_alignment/gene_alignment.hpp"
#include "alignment_reader/alignment_reader.hpp"
#include "statistics_estimator/statistics_estimator.hpp"
#include "statistics_exporter/statistics_exporter.hpp"

int shm_kmer_model_estimator::SHMkmerModelEstimator::Run() {
    const std::string boarder("=============");
    INFO(boarder << " SHM k-mer Model Calculator starts " << boarder);

    INFO(boarder << " Reading alignments starts " << boarder);
    INFO(std::string("input filename: ") << io_params_.input.input_filename);
    INFO(std::string("Strategy for checking: ") <<
        alignment_checker_params_.alignment_checker_method_names[
            static_cast<size_t> (alignment_checker_params_.alignment_checker_method)]);
    INFO(std::string("Strategy for cropping: ") <<
        alignment_cropper_params_.alignment_cropper_method_names[
            static_cast<size_t> (alignment_cropper_params_.alignment_cropper_method)]);

    ns_alignment_reader::AlignmentReader alignment_reader(io_params_.input.input_filename,
                                                          alignment_checker_params_,
                                                          alignment_cropper_params_);
    ns_gene_alignment::VectorReadGermlinePairs alignments(alignment_reader.read_alignments());
    INFO(boarder << " Reading alignments finishes " << boarder);

    INFO(boarder << " Estimating statistics starts " << boarder);
    INFO(std::string("Strategy for mutations: ") <<
         mutations_strategy_params_.mutation_strategy_method_names[
            static_cast<size_t> (mutations_strategy_params_.mutations_strategy_method)]);

    StatisticsEstimator statistics_estimator(mutations_strategy_params_);
    MutationsStatistics statistics =
        statistics_estimator.calculate_mutation_statistics(alignments);
    INFO(boarder << " Estimating statistics finishes" << boarder);

    INFO(boarder << " Exporting statistics starts " << boarder);
    INFO(std::string("output filename: ") << io_params_.output.output_filename);
    StatisticsExporter statistics_exporter(io_params_.output.output_filename);
    statistics_exporter.export_statistics(statistics);
    INFO(boarder << " Exporting statistics finishes " << boarder);

    return 0;
}