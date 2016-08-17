#pragma once

#include "single_loop_labelers/single_loop_helper.hpp"

namespace cdr_labeler {
    class BaseImmuneGeneCDRLabeler {
    public:
        virtual annotation_utils::CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene&) {
            return annotation_utils::CDRLabeling();
        }

        ~BaseImmuneGeneCDRLabeler() { }
    };

    typedef std::shared_ptr<BaseImmuneGeneCDRLabeler> BaseImmuneGeneCDRLabelerPtr;

    class DeNovoImmuneGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        germline_utils::ImmuneGeneType gene_type_;
        SingleLoopLabelerPtr cdr1_labeler_ptr_;
        SingleLoopLabelerPtr cdr2_labeler_ptr_;
        SingleLoopLabelerPtr cdr3_labeler_ptr_;

    public:
        DeNovoImmuneGeneCDRLabeler(germline_utils::ImmuneGeneType gene_type,
                                   const CDRLabelerConfig::CDRsParams &cdrs_params) :
                gene_type_(gene_type),
                cdr1_labeler_ptr_(SingleLoopLabelingHelper::GetCDR1Labeler(gene_type, cdrs_params)),
                cdr2_labeler_ptr_(SingleLoopLabelingHelper::GetCDR2Labeler(gene_type, cdrs_params)),
                cdr3_labeler_ptr_(SingleLoopLabelingHelper::GetCDR3Labeler(gene_type, cdrs_params)) { }

        annotation_utils::CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };
}