#include <verify.hpp>

#include "cdr_primitives.hpp"

namespace cdr_labeler {
    void DbCDRLabeling::AddGeneLabeling(const germline_utils::ImmuneGene &immune_gene, CDRLabeling labeling)  {
        cdr_labelings_.push_back(labeling);
        gene_name_index_map_[std::string(seqan::toCString(immune_gene.name()))] = cdr_labelings_.size() - 1;
        if(labeling.Empty())
            num_empty_labelings_++;
    }

    CDRLabeling DbCDRLabeling::GetLabelingByGene(const germline_utils::ImmuneGene &immune_gene) const {
        std::string gene_name = std::string(seqan::toCString(immune_gene.name()));
        VERIFY_MSG(gene_name_index_map_.find(gene_name) != gene_name_index_map_.end(),
                   "DB labeling does not contain information about gene " << immune_gene.name());
        return cdr_labelings_[gene_name_index_map_.at(gene_name)];
    }

    bool DbCDRLabeling::CDRLabelingIsEmpty(const germline_utils::ImmuneGene &immune_gene) const {
        auto labeling = GetLabelingByGene(immune_gene);
        return labeling.Empty();
    }
}