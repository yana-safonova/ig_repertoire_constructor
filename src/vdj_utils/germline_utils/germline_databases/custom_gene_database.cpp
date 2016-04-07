#include <verify.hpp>
#include <logger/logger.hpp>

#include "custom_gene_database.hpp"

namespace germline_utils {

    void CustomGeneDatabase::CheckConsistency(ImmuneGeneType gene_type) {
        VERIFY_MSG(gene_type.Segment() == segment_type_,
                   "Segment type " << gene_type.Segment() << " of new gene does not match with DB segment type: " <<
                   segment_type_);
    }

    void CustomGeneDatabase::AddDatabase(ImmuneGeneType gene_type, std::string filename) {
        CheckConsistency(gene_type);
        size_t db_index;
        if (gene_type_index_map_.find(gene_type) == gene_type_index_map_.end()) {
            gene_databases_.push_back(ImmuneGeneDatabase(gene_type));
            db_index = gene_databases_.size() - 1;
            gene_type_index_map_[gene_type] = db_index;
            immune_gene_types_.push_back(gene_type);
        }
        else
            db_index = gene_type_index_map_[gene_type];
        gene_databases_[db_index].AddGenesFromFile(filename);
        num_records_ += gene_databases_[db_index].size();
    }

    bool CustomGeneDatabase::ContainsImmuneGeneType(ImmuneGeneType gene_type) const {
        return gene_type_index_map_.find(gene_type) != gene_type_index_map_.end();
    }

    const ImmuneGeneDatabase& CustomGeneDatabase::GetDbByGeneType(ImmuneGeneType gene_type) const {
        VERIFY_MSG(ContainsImmuneGeneType(gene_type), "Custom DB does not contains gene type " << gene_type);
        return gene_databases_[gene_type_index_map_.at(gene_type)];
    }

    const ImmuneGene& CustomGeneDatabase::operator[](size_t index) const {
        return ImmuneGene();
    }
}