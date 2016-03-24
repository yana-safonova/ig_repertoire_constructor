#pragma once

#include "immune_gene_database.hpp"

namespace germline_utils {

// database stores databases for the specific type of gene segments (e.g., variable) and
// multiple lymphocyte and chain types
    class CustomGeneDatabase {
        SegmentType segment_type_;
        std::vector<ImmuneGeneDatabase> gene_databases_;
        // map from gene type from index in gene_databases vector
        std::unordered_map <ImmuneGeneType, size_t, ImmuneGeneTypeHasher> gene_type_index_map_;
        size_t num_records_;

        void CheckConsistency(ImmuneGeneType gene_type);

    public:
        CustomGeneDatabase(SegmentType segment_type) :
                segment_type_(segment_type),
                num_records_() { }

        void AddDatabase(ImmuneGeneType gene_type, std::string filename);

        size_t size() const { return num_records_; }

        SegmentType Segment() const { return segment_type_; }
    };
}