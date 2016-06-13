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
        std::vector<ImmuneGeneType> immune_gene_types_;
        size_t num_records_;
        // map from global index into a pair <db_index, record_local_index>
        std::unordered_map <size_t, std::pair<size_t, size_t>> gene_index_map_;

        void CheckConsistency(ImmuneGeneType gene_type);

        void UpdateGeneIndexMap(size_t db_index, size_t num_added_records);

        size_t AddImmuneGeneType(ImmuneGeneType gene_type);

    public:
        CustomGeneDatabase(SegmentType segment_type) :
                segment_type_(segment_type),
                num_records_() { }

        void AddDatabase(ImmuneGeneType gene_type, std::string filename);

        void AddImmuneGene(ImmuneGene immune_gene);

        size_t size() const { return num_records_; }

        size_t num_dbs() const { return gene_databases_.size(); }

        SegmentType Segment() const { return segment_type_; }

        typedef std::vector<ImmuneGeneType>::const_iterator ImmuneGeneTypeIter;

        ImmuneGeneTypeIter cbegin() const { return immune_gene_types_.cbegin(); }

        ImmuneGeneTypeIter cend() const { return immune_gene_types_.cend(); }

        bool ContainsImmuneGeneType(ImmuneGeneType gene_type) const;

        const ImmuneGeneDatabase& GetConstDbByGeneType(ImmuneGeneType gene_type) const;

        ImmuneGeneDatabase& GetDbByGeneType(ImmuneGeneType gene_type);

        const ImmuneGene& operator[](size_t index) const;
    };
}