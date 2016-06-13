#pragma once

#include <vj_alignment_info.hpp>
#include <annotation_utils/annotated_clone_set.hpp>

#include "germline_db_labeling.hpp"
#include "immune_gene_alignment_converter.hpp"


namespace cdr_labeler {
    class ReadCDRLabeler {
        const DbCDRLabeling& v_labeling_;
        const DbCDRLabeling& j_labeling_;

        vj_finder::ImmuneGeneAlignmentConverter alignment_converter_;

    public:
        ReadCDRLabeler(const DbCDRLabeling& v_labeling, const DbCDRLabeling& j_labeling) :
                v_labeling_(v_labeling),
                j_labeling_(j_labeling) { }

        annotation_utils::AnnotatedClone CreateAnnotatedClone(const vj_finder::VJHits &vj_hits);

        annotation_utils::CDRAnnotatedCloneSet CreateAnnotatedCloneSet(const vj_finder::VJAlignmentInfo &alignment_info);
    };
}