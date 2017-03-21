#pragma once

#include <germline_utils/germline_databases/immune_gene_database.hpp>
#include <annotation_utils/cdr_labeling_primitives.hpp>

namespace cdr_labeler {
    class SingleLoopLabeler {
    protected:
        germline_utils::ImmuneGeneType gene_type_;

        void CheckConsistencyFatal(const germline_utils::ImmuneGene &immune_gene);

        virtual annotation_utils::CDRRange ComputeRange(const germline_utils::ImmuneGene&,
                                                        annotation_utils::CDRRange) {
            return annotation_utils::CDRRange();
        }

    public:
        SingleLoopLabeler(germline_utils::ImmuneGeneType gene_type) : gene_type_(gene_type) { }

        annotation_utils::CDRRange ComputeLoopRange(const germline_utils::ImmuneGene &immune_gene,
                                                    annotation_utils::CDRRange previous_cdr);
    };

    typedef std::shared_ptr<SingleLoopLabeler> SingleLoopLabelerPtr;
}