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

    germline_utils::CustomGeneDatabase DbCDRLabeling::CreateFilteredDb() {
        germline_utils::CustomGeneDatabase filtered_db(germline_db_.Segment());
        for(auto it = germline_db_.cbegin(); it != germline_db_.cend(); it++) {
            auto specific_gene_db = germline_db_.GetDbByGeneType(*it);
            for(size_t i = 0; i < specific_gene_db.size(); i++) {
                auto gene_labeling = GetLabelingByGene(specific_gene_db[i]);
                if(!gene_labeling.Empty())
                    filtered_db.AddImmuneGene(specific_gene_db[i]);
            }
        }
        return filtered_db;
    }
}