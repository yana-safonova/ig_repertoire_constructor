#pragma once

#include "../../cdr_config.hpp"
#include "single_loop_labeler.hpp"

namespace cdr_labeler {
    class HCDR1Labeler : public SingleLoopLabeler {
        const CDRLabelerConfig::CDRsParams::HCDR1Params &params_;

        size_t ComputeStartPosition(const germline_utils::ImmuneGene &immune_gene);

        size_t ComputeEndPosition(const germline_utils::ImmuneGene &immune_gene, size_t start_pos);

        annotation_utils::CDRRange ComputeRange(const germline_utils::ImmuneGene &immune_gene,
                                                annotation_utils::CDRRange);

    public:
        HCDR1Labeler(const CDRLabelerConfig::CDRsParams::HCDR1Params &params) :
                SingleLoopLabeler(germline_utils::ImmuneGeneType(germline_utils::ChainType("IGH"),
                                                                 germline_utils::SegmentType::VariableSegment)),
                params_(params) { }
    };
}