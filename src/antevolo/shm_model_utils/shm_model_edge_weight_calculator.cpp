//
// Created by Andrew Bzikadze on 3/11/17.
//


#include <cmath>

#include "seqan/basic.h"
#include "shm_model_edge_weight_calculator.hpp"

using namespace shm_kmer_matrix_estimator;
using namespace annotation_utils;

namespace antevolo {

EvolutionaryEdgeAlignment
ShmModelEdgeWeightCalculator::get_prepared_strings(const BaseEvolutionaryEdge &edge) const
{
    auto &source_seqan = seqan::row(edge.SrcClone()->VAlignment().Alignment(), 1);
    auto &destination_seqan = seqan::row(edge.DstClone()->VAlignment().Alignment(), 1);

    std::string source(core::seqan_string_to_string(source_seqan));
    std::string destination(core::seqan_string_to_string(destination_seqan));

    auto source_it(source.begin());
    auto destination_it(destination.begin());
    while (source_it != source.end() and destination_it != destination.end()) {
        if (*source_it == '-' and *destination_it != '-') {
            destination.insert(destination_it, '-');
        } else if (*source_it != '-' and *destination_it == '-') {
            source.insert(source_it, '-');
        }
        ++source_it, ++destination_it;
    }

    size_t alignment_length(std::min(source.size(), destination.size()));
    source = source.substr(0, alignment_length);
    destination = destination.substr(0, alignment_length);
    std::string gene_id_src(core::seqan_string_to_string(edge.SrcClone()->VAlignment().subject().name()));
    std::string gene_id_dst(core::seqan_string_to_string(edge.DstClone()->VAlignment().subject().name()));

    VERIFY_MSG(gene_id_src == gene_id_dst,
               std::string("Gene id's should be the same. Src gene id = ") + gene_id_src +
               " , Dst gene is = " + gene_id_dst);

    size_t cdr1_start_src(edge.SrcClone()->CDR1Range().start_pos);
    size_t cdr1_end_src(edge.SrcClone()->CDR1Range().end_pos);
    size_t cdr2_start_src(edge.SrcClone()->CDR2Range().start_pos);
    size_t cdr2_end_src(edge.SrcClone()->CDR2Range().end_pos);

    size_t cdr1_start_dst(edge.DstClone()->CDR1Range().start_pos);
    size_t cdr1_end_dst(edge.DstClone()->CDR1Range().end_pos);
    size_t cdr2_start_dst(edge.DstClone()->CDR2Range().start_pos);
    size_t cdr2_end_dst(edge.DstClone()->CDR2Range().end_pos);

    VERIFY_MSG(cdr1_start_src == cdr1_start_dst,
               std::string("Cdr1 start should be equal. Src = ") + std::to_string(cdr1_start_src) +
               " , Dst = " + std::to_string(cdr1_start_dst));
    VERIFY_MSG(cdr2_start_src == cdr2_start_dst,
               std::string("Cdr2 start should be equal. Src = ") + std::to_string(cdr2_start_src) +
                   " , Dst = " + std::to_string(cdr2_start_dst));
    VERIFY_MSG(cdr1_end_src == cdr1_end_dst,
               std::string("Cdr1 end should be equal. Src = ") + std::to_string(cdr1_end_src) +
                   " , Dst = " + std::to_string(cdr1_end_dst));
    VERIFY_MSG(cdr2_end_src == cdr2_end_dst,
               std::string("Cdr2 end should be equal. Src = ") + std::to_string(cdr2_end_src) +
                   " , Dst = " + std::to_string(cdr2_end_dst));
    return EvolutionaryEdgeAlignment(std::move(source),
                                     std::move(destination),
                                     gene_id_src,
                                     edge.SrcClone()->HasStopCodon(),
                                     edge.SrcClone()->InFrame(),
                                     edge.SrcClone()->Productive(),
                                     cdr1_start_src,
                                     cdr1_end_src,
                                     cdr2_start_src,
                                     cdr2_end_src);
}

double ShmModelEdgeWeightCalculator::calculate_weigth_edge_per_position(const EvolutionaryEdgeAlignment &src_dst_pair,
                                                                        const size_t center_nucl_pos,
                                                                        const size_t kmer_len_) const
{
    std::string gene_substring =
        src_dst_pair.parent().substr(center_nucl_pos - kmer_len_ / 2, kmer_len_);
    const char center_nucl = src_dst_pair.son()[center_nucl_pos];

    if ((center_nucl == 'N') or
        (gene_substring.find_first_of('N') != std::string::npos))
        return 0;

    if ((src_dst_pair.son().find_first_of('-') != std::string::npos) or
        (src_dst_pair.parent().find_first_of('-') != std::string::npos))
        return 0;

    if ((center_nucl_pos >= src_dst_pair.cdr1_start() and center_nucl_pos <= src_dst_pair.cdr1_end()) or
        (center_nucl_pos >= src_dst_pair.cdr2_start() and center_nucl_pos <= src_dst_pair.cdr2_end())) {
        return model_.loglikelihood_kmer_nucl(gene_substring, center_nucl, StructuralRegion::CDR);
    }
    return model_.loglikelihood_kmer_nucl(gene_substring, center_nucl, StructuralRegion::FR);
}

double ShmModelEdgeWeightCalculator::calculate_weigth_edge(const BaseEvolutionaryEdge &edge) const {
    std::cout << "start calculating weight of " << edge.SrcNum() << " -> " << edge.DstNum() << std::endl;
    auto src_dst_pair(get_prepared_strings(edge));

    size_t kmer_len_(model_.kmer_len());
    std::vector<size_t>
        relevant_positions(mutation_strategy_->calculate_relevant_positions(src_dst_pair));

    double log_likelihood = 0;
    for (const auto& center_nucl_pos : relevant_positions) {
        double add_llklh = calculate_weigth_edge_per_position(src_dst_pair, center_nucl_pos, kmer_len_);
        if (not std::isnan(add_llklh)) {
            log_likelihood += add_llklh;
        }
    }
    std::cout << "end calculating weight of" << std::endl;
    return log_likelihood;
}

} // End namespace antevolo
