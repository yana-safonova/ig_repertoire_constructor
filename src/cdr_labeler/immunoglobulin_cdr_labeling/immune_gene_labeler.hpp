#pragma once

#include "hcdr1_labeler.hpp"
#include "hcdr2_labeler.hpp"
#include "hcdr3_v_labeler.hpp"
#include "hcdr3_j_labeler.hpp"

namespace cdr_labeler {
    class ImmuneGeneCDRHelper {
    public:
        static SingleLoopLabelerPtr GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);
    };

    class ImmuneGeneCDRLabeler {
        germline_utils::ImmuneGeneType gene_type_;
        SingleLoopLabelerPtr cdr1_labeler_ptr_;
        SingleLoopLabelerPtr cdr2_labeler_ptr_;
        SingleLoopLabelerPtr cdr3_labeler_ptr_;

    public:
        ImmuneGeneCDRLabeler(germline_utils::ImmuneGeneType gene_type,
                             const CDRLabelerConfig::CDRsParams &cdrs_params) :
                gene_type_(gene_type),
                cdr1_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR1Labeler(gene_type, cdrs_params)),
                cdr2_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR2Labeler(gene_type, cdrs_params)),
                cdr3_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR3Labeler(gene_type, cdrs_params)) { }

        CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };
}