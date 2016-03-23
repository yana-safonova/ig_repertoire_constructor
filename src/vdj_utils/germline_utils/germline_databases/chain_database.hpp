#pragma once

#include "immune_gene_database.hpp"

namespace germline_utils {

// database stores genes for specific type of chain, e.g., IGH or TRA
    class ChainDatabase {
        ChainType chain_type_;

        ImmuneGeneDatabase variable_genes_;
        ImmuneGeneDatabase diversity_genes_;
        ImmuneGeneDatabase join_genes_;

    public:
        ChainDatabase(ChainType chain_type) :
                chain_type_(chain_type),
                variable_genes_(ImmuneGeneType(chain_type, SegmentType::VariableSegment)),
                diversity_genes_(ImmuneGeneType(chain_type, SegmentType::DiversitySegment)),
                join_genes_(ImmuneGeneType(chain_type, SegmentType::JoinSegment)) { }

        void AddGenesFromFile(SegmentType segment_type, std::string filename);

        size_t GenesNumber(SegmentType segment_type) const;

        const ImmuneGene &GetGeneByIndex(SegmentType segment_type, size_t index) const;

        const ImmuneGeneDatabase &GetDb(SegmentType segment_type) const;

        ChainType Chain() const { return chain_type_; }
    };

    std::ostream &operator<<(std::ostream &out, const ChainDatabase &chain_db);
}
