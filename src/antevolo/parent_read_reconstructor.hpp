#pragma once

#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {

    class ParentReadReconstructor {
    public:
        static std::tuple<core::Read,
                seqan::Align<seqan::Dna5String>,
                seqan::Align<seqan::Dna5String>>
        ReconstructParentRead(
               const annotation_utils::AnnotatedClone& clone1,
               const annotation_utils::AnnotatedClone& clone2,
               size_t id,
               size_t gene_cdr3_start_pos,
               size_t gene_cdr3_end_pos);

    private:
        using TRow =  seqan::Row<seqan::Align<seqan::Dna5String, seqan::ArrayGaps>>::Type;

        static std::vector<size_t> TraverseAlignments(
                const TRow& alignment1,
                const TRow& alignment2,
                const TRow& gene_alignment, //check for equality
                size_t read1_start_pos,
                size_t read2_start_pos,
                size_t read1_end_pos,
                size_t read2_end_pos,
                std::string &res_string);
    };
}
