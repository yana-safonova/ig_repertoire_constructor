//
// Created by Andrew Bzikadze on 5/23/16.
//

#include <string>
#include <seqan/file.h>

#include "trivial_strategy.hpp"

using namespace ns_abstract_mutation_strategy;
using namespace ns_gene_alignment;


std::vector<size_t> TrivialMutationStrategy::calculate_relevant_positions
    (ns_gene_alignment::ReadGermlineAlignment &alignment) const {
    std::vector<size_t> relevant_positions(alignment.read().size() - kmer_len_ + 1);
    std::iota(relevant_positions.begin(), relevant_positions.end(), kmer_len_ / 2);
    return relevant_positions;
}

