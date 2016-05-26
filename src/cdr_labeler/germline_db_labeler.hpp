#pragma once

#include "cdr_primitives.hpp"

namespace cdr_labeler {
    class GermlineDbLabeler {
    protected:
        const germline_utils::CustomGeneDatabase &gene_db_;

    public:
        GermlineDbLabeler(const germline_utils::CustomGeneDatabase &gene_db) :
                gene_db_(gene_db) { }

        DbCDRLabeling ComputeLabeling();
    };
}