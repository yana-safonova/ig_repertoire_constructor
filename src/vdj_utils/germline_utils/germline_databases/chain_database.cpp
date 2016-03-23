#include <verify.hpp>

#include "chain_database.hpp"

namespace germline_utils {

// todo: prohibit addition of diversity genes for VJ chains
    void ChainDatabase::AddGenesFromFile(SegmentType segment_type, std::string filename) {
        if (segment_type == SegmentType::VariableSegment)
            variable_genes_.AddGenesFromFile(filename);
        else if (segment_type == SegmentType::DiversitySegment)
            diversity_genes_.AddGenesFromFile(filename);
        else if (segment_type == SegmentType::JoinSegment)
            join_genes_.AddGenesFromFile(filename);
    }

    size_t ChainDatabase::GenesNumber(SegmentType segment_type) const {
        if (segment_type == SegmentType::VariableSegment)
            return variable_genes_.size();
        if (segment_type == SegmentType::DiversitySegment)
            return diversity_genes_.size();
        if (segment_type == SegmentType::JoinSegment)
            return join_genes_.size();
        return 0;
    }

    const ImmuneGene &ChainDatabase::GetGeneByIndex(SegmentType segment_type, size_t index) const {
        VERIFY_MSG(segment_type != SegmentType::UnknownSegment, "Requested segment type is unknown");
        if (segment_type == SegmentType::VariableSegment)
            return variable_genes_[index];
        if (segment_type == SegmentType::DiversitySegment)
            return diversity_genes_[index];
        return join_genes_[index];
    }

    const ImmuneGeneDatabase &ChainDatabase::GetDb(SegmentType segment_type) const {
        VERIFY_MSG(segment_type != SegmentType::UnknownSegment, "Requested segment type is unknown");
        if (segment_type == SegmentType::VariableSegment)
            return variable_genes_;
        if (segment_type == SegmentType::DiversitySegment)
            return diversity_genes_;
        return join_genes_;
    }

    std::ostream &operator<<(std::ostream &out, const ChainDatabase &chain_db) {
        out << "Database for immune chain " << chain_db.Chain() << std::endl;
        if (chain_db.GenesNumber(SegmentType::VariableSegment) != 0)
            out << chain_db.GetDb(SegmentType::VariableSegment) << std::endl;
        if (chain_db.GenesNumber(SegmentType::DiversitySegment) != 0)
            out << chain_db.GetDb(SegmentType::DiversitySegment) << std::endl;
        if (chain_db.GenesNumber(SegmentType::JoinSegment) != 0)
            out << chain_db.GetDb(SegmentType::JoinSegment);
        return out;
    }
}
