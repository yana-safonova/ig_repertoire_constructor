#pragma once

#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {

    class IntersectionParentConstructor {
    public:
        static annotation_utils::AnnotatedClone ReconstructParent(
                const std::shared_ptr<annotation_utils::AnnotatedClone> &clone1,
                const std::shared_ptr<annotation_utils::AnnotatedClone> &clone2);

    private:
        static void TraverseReads(const seqan::String& read1,
                                  const seqan::String& read2,
                                  size_t read1_last_pos,
                                  size_t read2_last_pos,
                                  size_t read1_end_pos,
                                  size_t read2_end_pos,
                                  annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it1,
                                  annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it2,
                                  annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end1,
                                  annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end2,
                                  std::string& res_string,
                                  annotation_utils::CDRRange& res_cdr1_range,
                                  annotation_utils::CDRRange& res_cdr2_range,
                                  annotation_utils::CDRRange& res_cdr3_range);

    };

}
