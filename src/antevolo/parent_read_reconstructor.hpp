#pragma once

#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {

    class ParentReadReconstructor {
    public:
        static core::Read ReconstructParentRead(
                const std::shared_ptr<annotation_utils::AnnotatedClone> &clone1,
                const std::shared_ptr<annotation_utils::AnnotatedClone> &clone2,
                size_t id);

    private:
        static void TraverseReads(const seqan::Dna5String& read1,
                           const seqan::Dna5String& read2,
                           size_t read1_last_pos,
                           size_t read2_last_pos,
                           size_t read1_end_pos,
                           size_t read2_end_pos,
                           annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it1,
                           annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it2,
                           annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end1,
                           annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end2,
                           std::string& res_string);
    };
}
