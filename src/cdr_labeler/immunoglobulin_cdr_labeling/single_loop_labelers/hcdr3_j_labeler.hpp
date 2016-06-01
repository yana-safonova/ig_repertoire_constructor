#pragma once

#include "../../cdr_config.hpp"
#include "single_loop_labeler.hpp"

namespace cdr_labeler {
    class HCDR3JLabeler : public SingleLoopLabeler {
        const CDRLabelerConfig::CDRsParams::HCDR3Params &params_;

        size_t ComputeEndPosition(const germline_utils::ImmuneGene &immune_gene);

        CDRRange ComputeRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr);

    public:
        HCDR3JLabeler(const CDRLabelerConfig::CDRsParams::HCDR3Params &params) :
                SingleLoopLabeler(germline_utils::ImmuneGeneType(germline_utils::ChainType("IGH"),
                                                                 germline_utils::SegmentType::JoinSegment)),
                params_(params) { }
    };
}