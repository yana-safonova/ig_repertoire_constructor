#include <verify.hpp>

#include "chain_database.hpp"

namespace germline_utils {
    void ChainDatabase::InitializeImmuneGeneDbs() {
        immune_gene_dbs_[SegmentType::VariableSegment] =
                ImmuneGeneDatabase(ImmuneGeneType(chain_type_, SegmentType::VariableSegment));
        immune_gene_dbs_[SegmentType::JoinSegment] =
                ImmuneGeneDatabase(ImmuneGeneType(chain_type_, SegmentType::JoinSegment));
        if(chain_type_.IsVDJ())
            immune_gene_dbs_[SegmentType::DiversitySegment] =
                    ImmuneGeneDatabase(ImmuneGeneType(chain_type_, SegmentType::DiversitySegment));
    }

    void ChainDatabase::CheckConsistency(SegmentType segment_type) const {
        if(segment_type == SegmentType::DiversitySegment)
        VERIFY_MSG(chain_type_.IsVDJ(), "Chain " << chain_type_ << " does not contain D genes");
        VERIFY_MSG(segment_type != SegmentType::UnknownSegment, "Requested segment type is unknown");
    }

    void ChainDatabase::AddGenesFromFile(SegmentType segment_type, std::string filename) {
        CheckConsistency(segment_type);
        immune_gene_dbs_.at(segment_type).AddGenesFromFile(filename);
    }

    size_t ChainDatabase::GenesNumber(SegmentType segment_type) const {
        CheckConsistency(segment_type);
        return immune_gene_dbs_.at(segment_type).size();
    }

    const ImmuneGene &ChainDatabase::GetGeneByIndex(SegmentType segment_type, size_t index) const {
        CheckConsistency(segment_type);
        return immune_gene_dbs_.at(segment_type)[index];
    }

    const ImmuneGeneDatabase &ChainDatabase::GetDb(SegmentType segment_type) const {
        CheckConsistency(segment_type);
        return immune_gene_dbs_.at(segment_type);
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
