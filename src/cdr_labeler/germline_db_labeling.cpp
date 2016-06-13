#include <verify.hpp>

#include "germline_db_labeling.hpp"

namespace cdr_labeler {
    using namespace annotation_utils;

    bool DbCDRLabeling::LabelingIsValid(const germline_utils::ImmuneGene &immune_gene,
                                        annotation_utils::CDRLabeling labeling) const {
        size_t gene_length = immune_gene.length();
        // todo: refactor it!
        if(immune_gene.GeneType().Segment() == germline_utils::SegmentType::VariableSegment)
            return labeling.cdr1.start_pos < gene_length and labeling.cdr1.end_pos < gene_length and
                    labeling.cdr2.start_pos < gene_length and labeling.cdr2.end_pos < gene_length and
                    labeling.cdr3.start_pos < gene_length;
        return labeling.cdr3.end_pos < gene_length;
    }

    void DbCDRLabeling::AddGeneLabeling(const germline_utils::ImmuneGene &immune_gene, CDRLabeling labeling)  {
        cdr_labelings_.push_back(labeling);
        gene_name_index_map_[std::string(seqan::toCString(immune_gene.name()))] = cdr_labelings_.size() - 1;
        if(labeling.Empty() or !LabelingIsValid(immune_gene, labeling))
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
            auto specific_gene_db = germline_db_.GetConstDbByGeneType(*it);
            for(size_t i = 0; i < specific_gene_db.size(); i++) {
                auto gene_labeling = GetLabelingByGene(specific_gene_db[i]);
                if(!gene_labeling.Empty() and LabelingIsValid(specific_gene_db[i], gene_labeling))
                    filtered_db.AddImmuneGene(specific_gene_db[i]);
            }
        }
        return filtered_db;
    }
}