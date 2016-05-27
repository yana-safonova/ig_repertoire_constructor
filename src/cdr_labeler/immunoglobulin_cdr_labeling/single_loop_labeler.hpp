#pragma once

#include "../cdr_primitives.hpp"

namespace cdr_labeler {
    class SingleLoopLabeler {
    protected:
        germline_utils::ImmuneGeneType gene_type_;

        void CheckConsistencyFatal(const germline_utils::ImmuneGene &immune_gene);

        virtual CDRRange ComputeRange(const germline_utils::ImmuneGene &immune_gene,
                                      CDRRange) {
            return CDRRange();
        }

    public:
        SingleLoopLabeler(germline_utils::ImmuneGeneType gene_type) : gene_type_(gene_type) { }

        CDRRange ComputeLoopRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr);
    };

    typedef std::shared_ptr<SingleLoopLabeler> SingleLoopLabelerPtr;
}