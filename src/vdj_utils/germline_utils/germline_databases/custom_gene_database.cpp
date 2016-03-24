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
        }
        else
            db_index = gene_type_index_map_[gene_type];
        gene_databases_[db_index].AddGenesFromFile(filename);
        num_records_ += gene_databases_[db_index].size();
    }
}