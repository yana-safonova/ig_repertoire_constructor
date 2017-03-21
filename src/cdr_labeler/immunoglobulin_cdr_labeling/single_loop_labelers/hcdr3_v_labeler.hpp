#pragma once

#include "../../cdr_config.hpp"
#include "single_loop_labeler.hpp"

namespace cdr_labeler {
    class HCDR3VLabeler : public SingleLoopLabeler {
        //const CDRLabelerConfig::CDRsParams::HCDR3Params &params_;

        size_t ComputeStartPosition(const germline_utils::ImmuneGene &immune_gene,
                                    annotation_utils::CDRRange previous_cdr);

        annotation_utils::CDRRange ComputeRange(const germline_utils::ImmuneGene &immune_gene,
                                                annotation_utils::CDRRange previous_cdr);

    public:
        HCDR3VLabeler(const CDRLabelerConfig::CDRsParams::HCDR3Params &) :
                SingleLoopLabeler(germline_utils::ImmuneGeneType(germline_utils::ChainType("IGH"),
                                                                 germline_utils::SegmentType::VariableSegment))
                /*, params_(params)*/ { }
    };
}