#pragma once

#include "hcdr1_labeler.hpp"
#include "hcdr2_labeler.hpp"
#include "hcdr3_v_labeler.hpp"
#include "hcdr3_j_labeler.hpp"

namespace cdr_labeler {
    class SingleLoopLabelingHelper {
    public:
        static SingleLoopLabelerPtr GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);
    };
}