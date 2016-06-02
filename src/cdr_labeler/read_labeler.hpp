#pragma once

#include "germline_db_labeling.hpp"
#include <annotation_utils/cdr_annotated_clone.hpp>
#include <vj_alignment_structs.hpp>

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

        annotation_utils::CDRAnnotatedClone CreateAnnotatedClone(const vj_finder::VJHits &vj_hits);
    };
}