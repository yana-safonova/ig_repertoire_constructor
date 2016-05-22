//
// Created by Andrew Bzikadze on 5/18/16.
//

#pragma once

#include <string>
#include <vector>

namespace ns_gene_alignment {
    // using DnaGapped = seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'-'>>;
    // using DnaGappedString = seqan::String<DnaGapped>;
    // using DnaGappedAlignment = seqan::Align<DnaGappedString, seqan::ArrayGaps>;

    // We place germline and read as pure c++ strings here for now.
    using GermlineReadPair = std::pair<std::string, std::string>;
    using VectorGermlineReadPairs = std::vector<ns_gene_alignment::GermlineReadPair>;
}
