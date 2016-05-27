#pragma once

#include "cdr_primitives.hpp"
#include "cdr_config.hpp"

namespace cdr_labeler {
    class GermlineDbLabeler {
    protected:
        const germline_utils::CustomGeneDatabase &gene_db_;
        const CDRLabelerConfig::CDRsParams &cdr_params_;

    public:
        GermlineDbLabeler(const germline_utils::CustomGeneDatabase &gene_db,
                          const CDRLabelerConfig::CDRsParams &cdr_params) :
                gene_db_(gene_db), cdr_params_(cdr_params) { }

        DbCDRLabeling ComputeLabeling();
    };
}