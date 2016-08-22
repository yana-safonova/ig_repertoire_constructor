#pragma once

#include <unordered_map>

#include "immune_gene_database.hpp"

namespace germline_utils {

// database stores genes for specific type of chain, e.g., IGH or TRA
    class ChainDatabase {
        ChainType chain_type_;

        std::unordered_map<SegmentType, ImmuneGeneDatabase, SegmentTypeHasher> immune_gene_dbs_;

        void InitializeImmuneGeneDbs();

        void CheckConsistency(SegmentType segment_type) const;

    public:
        ChainDatabase(ChainType chain_type) :
                chain_type_(chain_type) {
            InitializeImmuneGeneDbs();
        }

        void AddGenesFromFile(SegmentType segment_type, std::string filename);

        size_t GenesNumber(SegmentType segment_type) const;

        const ImmuneGene &GetGeneByIndex(SegmentType segment_type, size_t index) const;

        const ImmuneGeneDatabase &GetDb(SegmentType segment_type) const;

        ChainType Chain() const { return chain_type_; }

        bool IsVDJ() const { return chain_type_.IsVDJ(); }

        bool ContainsGeneSegments(SegmentType segment_type) const;

        bool Complete() const;
    };

    std::ostream &operator<<(std::ostream &out, const ChainDatabase &chain_db);
}
