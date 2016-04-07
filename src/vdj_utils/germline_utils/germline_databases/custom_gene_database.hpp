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

        void CheckConsistency(ImmuneGeneType gene_type);

    public:
        CustomGeneDatabase(SegmentType segment_type) :
                segment_type_(segment_type),
                num_records_() { }

        void AddDatabase(ImmuneGeneType gene_type, std::string filename);

        size_t size() const { return num_records_; }

        SegmentType Segment() const { return segment_type_; }

        // todo: implement me and refactor
        size_t GetSequenceLengthByIndex(size_t index) const { return 0; }

        typedef std::vector<ImmuneGeneType>::const_iterator ImmuneGeneTypeIter;

        ImmuneGeneTypeIter cbegin() const { return immune_gene_types_.cbegin(); }

        ImmuneGeneTypeIter cend() const { return immune_gene_types_.cend(); }

        bool ContainsImmuneGeneType(ImmuneGeneType gene_type) const;

        const ImmuneGeneDatabase& GetDbByGeneType(ImmuneGeneType gene_type) const;

        const ImmuneGene& operator[](size_t index) const;
    };
}