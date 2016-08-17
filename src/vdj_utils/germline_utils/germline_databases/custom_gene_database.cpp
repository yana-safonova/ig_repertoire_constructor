#include <verify.hpp>
#include <logger/logger.hpp>

#include "custom_gene_database.hpp"

namespace germline_utils {

    void CustomGeneDatabase::CheckConsistency(ImmuneGeneType gene_type) {
        VERIFY_MSG(gene_type.Segment() == segment_type_,
                   "Segment type " << gene_type.Segment() << " of new gene does not match with DB segment type: " <<
                   segment_type_);
    }

    void CustomGeneDatabase::UpdateGeneIndexMap(size_t db_index, size_t num_added_records) {
        size_t num_db_records = gene_databases_[db_index].size();
        for(size_t i = 0; i < num_added_records; i++) {
            size_t global_id = num_records_ - num_added_records + i;
            size_t local_id = num_db_records - num_added_records + i;
            gene_index_map_[global_id] = std::make_pair(db_index, local_id);
        }
    }

    size_t CustomGeneDatabase::AddImmuneGeneType(ImmuneGeneType gene_type) {
        gene_databases_.push_back(ImmuneGeneDatabase(gene_type));
        size_t db_index = gene_databases_.size() - 1;
        gene_type_index_map_[gene_type] = db_index;
        immune_gene_types_.push_back(gene_type);
        return db_index;
    }

    void CustomGeneDatabase::AddDatabase(ImmuneGeneType gene_type, std::string filename) {
        CheckConsistency(gene_type);
        size_t db_index;
        if (!ContainsImmuneGeneType(gene_type))
            db_index = AddImmuneGeneType(gene_type);
        else
            db_index = gene_type_index_map_[gene_type];
        size_t num_added_records = gene_databases_[db_index].AddGenesFromFile(filename);
        num_records_ += num_added_records;
        UpdateGeneIndexMap(db_index, num_added_records);
    }

    void CustomGeneDatabase::AddImmuneGene(ImmuneGene immune_gene) {
        CheckConsistency(immune_gene.GeneType());
        size_t db_index;
        if (!ContainsImmuneGeneType(immune_gene.GeneType()))
            db_index = AddImmuneGeneType(immune_gene.GeneType());
        else
            db_index = gene_type_index_map_[immune_gene.GeneType()];
        gene_databases_[db_index].AddImmuneGene(immune_gene);
        num_records_++;
        UpdateGeneIndexMap(db_index, 1);
    }

    bool CustomGeneDatabase::ContainsImmuneGeneType(ImmuneGeneType gene_type) const {
        return gene_type_index_map_.find(gene_type) != gene_type_index_map_.end();
    }

    const ImmuneGeneDatabase& CustomGeneDatabase::GetConstDbByGeneType(ImmuneGeneType gene_type) const {
        VERIFY_MSG(ContainsImmuneGeneType(gene_type), "Custom DB does not contains gene type " << gene_type);
        return gene_databases_.at(gene_type_index_map_.at(gene_type));
    }

    ImmuneGeneDatabase& CustomGeneDatabase::GetDbByGeneType(ImmuneGeneType gene_type) {
        VERIFY_MSG(ContainsImmuneGeneType(gene_type), "Custom DB does not contains gene type " << gene_type);
        return gene_databases_.at(gene_type_index_map_.at(gene_type));
    }

    const ImmuneGene& CustomGeneDatabase::operator[](size_t index) const {
        VERIFY_MSG(gene_index_map_.find(index) != gene_index_map_.end(), "Index " << index <<
                " is not presented in gene map");
        auto index_pair = gene_index_map_.at(index);
        return gene_databases_[index_pair.first][index_pair.second];
    }
}