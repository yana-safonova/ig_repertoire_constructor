#include "read_labeler.hpp"
#include <annotation_utils/annotated_clone_calculator.hpp>

namespace cdr_labeler {
    std::shared_ptr<annotation_utils::BaseAACalculator> ReadCDRLabeler::GetAACalculator() {
        return std::shared_ptr<annotation_utils::BaseAACalculator>(new annotation_utils::SimpleAACalculator());
    }

    std::shared_ptr<annotation_utils::BaseSHMCalculator> ReadCDRLabeler::GetVSHMCalculator() {
        if(shm_config_.shm_finding_algorithm ==
                CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::FilteringSHMAlgorithm)
            return std::shared_ptr<annotation_utils::BaseSHMCalculator>(
                    new annotation_utils::StartEndFilteringSHMCalculator(shm_config_.shm_filtering_params.v_start_max_skipped,
                                                                         shm_config_.shm_filtering_params.v_end_max_skipped));
        if(shm_config_.shm_finding_algorithm ==
                    CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm)
            return std::shared_ptr<annotation_utils::BaseSHMCalculator>(
                    new annotation_utils::CDRFilteringSHMCalculator());
        VERIFY_MSG(false, "Unknown SHM finding algorithm");
        return std::shared_ptr<annotation_utils::BaseSHMCalculator>(NULL);
    }

    std::shared_ptr<annotation_utils::BaseSHMCalculator> ReadCDRLabeler::GetJSHMCalculator() {
        if(shm_config_.shm_finding_algorithm ==
           CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::FilteringSHMAlgorithm)
            return std::shared_ptr<annotation_utils::BaseSHMCalculator>(
                    new annotation_utils::StartEndFilteringSHMCalculator(shm_config_.shm_filtering_params.j_start_max_skipped,
                                                                         shm_config_.shm_filtering_params.j_end_max_skipped));
        if(shm_config_.shm_finding_algorithm ==
           CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm)
            return std::shared_ptr<annotation_utils::BaseSHMCalculator>(
                    new annotation_utils::CDRFilteringSHMCalculator());
        VERIFY_MSG(false, "Unknown SHM finding algorithm");
        return std::shared_ptr<annotation_utils::BaseSHMCalculator>(NULL);
    }

    annotation_utils::AnnotatedClone ReadCDRLabeler::CreateAnnotatedClone(const vj_finder::VJHits &vj_hits) {
        auto v_hit = vj_hits.GetVHitByIndex(0);
        auto v_alignment = alignment_converter_.ConvertToAlignment(v_hit.ImmuneGene(),
                                                                   vj_hits.Read(),
                                                                   v_hit.BlockAlignment());
        auto v_cdr_labeling = v_labeling_.GetLabelingByGene(v_hit.ImmuneGene());
        annotation_utils::CDRRange read_cdr1(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr1.start_pos),
                                             v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr1.end_pos));
        annotation_utils::CDRRange read_cdr2(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr2.start_pos),
                                             v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr2.end_pos));
        auto j_hit = vj_hits.GetJHitByIndex(0);
        auto j_alignment = alignment_converter_.ConvertToAlignment(j_hit.ImmuneGene(),
                                                                   vj_hits.Read(),
                                                                   j_hit.BlockAlignment());
        auto j_cdr_labeling = j_labeling_.GetLabelingByGene(j_hit.ImmuneGene());
        annotation_utils::CDRRange read_cdr3(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr3.start_pos),
                                             j_alignment.QueryPositionBySubjectPosition(j_cdr_labeling.cdr3.end_pos));

        return clone_calculator_.ComputeAnnotatedClone(vj_hits.Read(),
                                                       annotation_utils::CDRLabeling(read_cdr1, read_cdr2, read_cdr3),
                                                       v_alignment, j_alignment);
    }

    annotation_utils::CDRAnnotatedCloneSet ReadCDRLabeler::CreateAnnotatedCloneSet(
            const vj_finder::VJAlignmentInfo &alignment_info) {
        annotation_utils::CDRAnnotatedCloneSet clone_set;
        for(size_t i = 0; i < alignment_info.NumVJHits(); i++) {
            auto vj_hit = alignment_info.GetVJHitsByIndex(i);
            clone_set.AddClone(CreateAnnotatedClone(vj_hit));
        }
        INFO(clone_set.size() << " annotated sequences were created");
        return clone_set;
    }
}