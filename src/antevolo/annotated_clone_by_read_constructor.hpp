#pragma once


#include <annotation_utils/annotated_clone.hpp>
#include <read_labeler.hpp>
#include <vj_finder_config.hpp>

namespace antevolo {

    class AnnotatedCloneByReadConstructor {
        const germline_utils::CustomGeneDatabase& labeled_v_db_;
        const germline_utils::CustomGeneDatabase& labeled_j_db_;
        const vj_finder::VJFinderConfig::AlgorithmParams& vj_finder_params_;
        cdr_labeler::ReadCDRLabeler& read_labeler_;
    public:
        AnnotatedCloneByReadConstructor(const germline_utils::CustomGeneDatabase& labeled_v_db,
                                      const germline_utils::CustomGeneDatabase& labeled_j_db,
                                      const vj_finder::VJFinderConfig::AlgorithmParams& vj_finder_params,
                                      cdr_labeler::ReadCDRLabeler& read_labeler) :
                labeled_v_db_(labeled_v_db),
                labeled_j_db_(labeled_j_db),
                vj_finder_params_(vj_finder_params),
                read_labeler_(read_labeler) {}

        annotation_utils::AnnotatedClone GetCloneByRead(core::Read& read) const;

    };

}
