#pragma once

#include "../../cdr_config.hpp"
#include "single_loop_labeler.hpp"

namespace cdr_labeler {
    class HCDR2Labeler : public SingleLoopLabeler {
        const CDRLabelerConfig::CDRsParams::HCDR2Params &params_;

        size_t ComputeStartPosition(const germline_utils::ImmuneGene &immune_gene,
                                    annotation_utils::CDRRange previous_cdr);

        size_t ComputeEndPosition(const germline_utils::ImmuneGene &immune_gene, size_t start_pos);

        annotation_utils::CDRRange ComputeRange(const germline_utils::ImmuneGene &immune_gene,
                                                annotation_utils::CDRRange previous_cdr);

    public:
        HCDR2Labeler(const CDRLabelerConfig::CDRsParams::HCDR2Params &params) :
                SingleLoopLabeler(germline_utils::ImmuneGeneType(germline_utils::ChainType("IGH"),
                                                                 germline_utils::SegmentType::VariableSegment)),
                params_(params) { }
    };
}