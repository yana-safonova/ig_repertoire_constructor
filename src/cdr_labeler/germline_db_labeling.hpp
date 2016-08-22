#pragma once

#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include <annotation_utils/cdr_labeling_primitives.hpp>

namespace cdr_labeler {
    class DbCDRLabeling {
        const germline_utils::CustomGeneDatabase &germline_db_;

        std::vector<annotation_utils::CDRLabeling> cdr_labelings_;
        std::unordered_map<std::string, size_t> gene_name_index_map_;

        size_t num_empty_labelings_;

        bool LabelingIsValid(const germline_utils::ImmuneGene &immune_gene,
                                annotation_utils::CDRLabeling labeling) const;

    public:
        DbCDRLabeling(const germline_utils::CustomGeneDatabase &germline_db) : germline_db_(germline_db),
                                                                               num_empty_labelings_() { }

        void AddGeneLabeling(const germline_utils::ImmuneGene &immune_gene,
                             annotation_utils::CDRLabeling labeling);

        annotation_utils::CDRLabeling GetLabelingByGene(const germline_utils::ImmuneGene &immune_gene) const;

        bool CDRLabelingIsEmpty(const germline_utils::ImmuneGene &immune_gene) const;

        size_t NumEmptyLabelings() const { return num_empty_labelings_; }

        germline_utils::CustomGeneDatabase CreateFilteredDb();
    };
}