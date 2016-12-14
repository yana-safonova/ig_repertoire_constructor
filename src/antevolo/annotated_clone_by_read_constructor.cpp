#include <vj_query_aligner.hpp>
#include "annotated_clone_by_read_constructor.hpp"

namespace antevolo {
    annotation_utils::AnnotatedClone AnnotatedCloneByReadConstructor::GetCloneByRead(core::Read& read) const {

        vj_finder::VJQueryAligner vj_query_aligner(vj_finder_params_, labeled_v_db_, labeled_j_db_);
        vj_finder::VJHits vj_hits = vj_query_aligner.Align(read);
        return read_labeler_.CreateAnnotatedClone(vj_hits);
    }
}
