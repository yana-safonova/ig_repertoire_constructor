#pragma once

#include <vj_alignment_info.hpp>
#include <annotation_utils/annotated_clone_set.hpp>
#include <annotation_utils/annotated_clone_calculator.hpp>

#include "germline_db_labeling.hpp"
#include "immune_gene_alignment_converter.hpp"
#include "cdr_config.hpp"


namespace cdr_labeler {
    class ReadCDRLabeler {
        const CDRLabelerConfig::SHMFindingParams &shm_config_;
        const DbCDRLabeling& v_labeling_;
        const DbCDRLabeling& j_labeling_;

        annotation_utils::AnnotatedCloneCalculator clone_calculator_;
        vj_finder::ImmuneGeneAlignmentConverter alignment_converter_;

        std::shared_ptr<annotation_utils::BaseAACalculator> GetAACalculator();

        std::shared_ptr<annotation_utils::BaseSHMCalculator> GetVSHMCalculator();

        std::shared_ptr<annotation_utils::BaseSHMCalculator> GetJSHMCalculator();

    public:
        ReadCDRLabeler(const CDRLabelerConfig::SHMFindingParams &shm_config,
                const DbCDRLabeling& v_labeling, const DbCDRLabeling& j_labeling) :
                shm_config_(shm_config),
                v_labeling_(v_labeling),
                j_labeling_(j_labeling),
                clone_calculator_(GetAACalculator(), GetVSHMCalculator(), GetJSHMCalculator()){ }

        annotation_utils::AnnotatedClone CreateAnnotatedClone(const vj_finder::VJHits &vj_hits);

        annotation_utils::CDRAnnotatedCloneSet CreateAnnotatedCloneSet(const vj_finder::VJAlignmentInfo &alignment_info);
    };
}