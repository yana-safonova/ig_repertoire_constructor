#pragma once

#include <block_alignment/block_alignment_converter.hpp>

namespace vj_finder {
    class ImmuneGeneAlignmentConverter : public algorithms::BlockAlignmentConverter<
            germline_utils::ImmuneGene, core::Read> {
    public:
        seqan::Dna5String GetSubjectString(const germline_utils::ImmuneGene &subject) {
            return subject.seq();
        }

        seqan::Dna5String GetQueryString(const core::Read& query) {
            return query.seq;
        }
    };
}