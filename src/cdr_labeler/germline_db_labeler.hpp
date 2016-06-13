#pragma once

#include "germline_db_labeling.hpp"
#include "cdr_config.hpp"

#include "immunoglobulin_cdr_labeling/immune_gene_labeler.hpp"

namespace cdr_labeler {
    class GermlineDbLabeler {
    protected:
        germline_utils::CustomGeneDatabase &gene_db_;
        const CDRLabelerConfig::CDRsParams &cdr_params_;

        std::shared_ptr<BaseImmuneGeneCDRLabeler> GetImmuneGeneLabeler(germline_utils::ImmuneGeneType gene_type);

    public:
        GermlineDbLabeler(germline_utils::CustomGeneDatabase &gene_db,
                          const CDRLabelerConfig::CDRsParams &cdr_params) :
                gene_db_(gene_db), cdr_params_(cdr_params) { }

        DbCDRLabeling ComputeLabeling();
    };
}