#pragma once

#include "hcdr1_labeler.hpp"
#include "hcdr2_labeler.hpp"

namespace cdr_labeler {
    class ImmuneGeneCDRHelper {
    public:
        static std::shared_ptr<SingleLoopLabeler> GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                 const CDRLabelerConfig::CDRsParams &cdr_params);

        static std::shared_ptr<SingleLoopLabeler> GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                 const CDRLabelerConfig::CDRsParams &cdr_params);

        //static std::shared_ptr<SingleLoopLabeler> GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
        //                                                         const CDRLabelerConfig::CDRsParams &cdr_params);
    };

    class ImmuneGeneCDRLabeler {
        germline_utils::ImmuneGeneType gene_type_;
        std::shared_ptr<SingleLoopLabeler> cdr1_labeler_ptr_;
        std::shared_ptr<SingleLoopLabeler> cdr2_labeler_ptr_;
        //std::shared_ptr<SingleLoopLabeler> cdr3_labeler_ptr_;

    public:
        ImmuneGeneCDRLabeler(germline_utils::ImmuneGeneType gene_type,
                             const CDRLabelerConfig::CDRsParams &cdrs_params) :
                gene_type_(gene_type),
                cdr1_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR1Labeler(gene_type, cdrs_params)),
                cdr2_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR2Labeler(gene_type, cdrs_params)) { }
        //        cdr3_labeler_(gene_type) { }

        CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };
}