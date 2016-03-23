#include <verify.hpp>

#include "germline_gene_type.hpp"

namespace germline_utils {

    std::ostream &operator<<(std::ostream &out, const SegmentType &segment_type) {
        if (segment_type == SegmentType::VariableSegment)
            out << "V";
        else if (segment_type == SegmentType::DiversitySegment)
            out << "D";
        else if (segment_type == SegmentType::JoinSegment)
            out << "J";
        else
            out << "Unknown gene segment";
        return out;
    }

    bool ImmuneGeneType::CheckChainSegmentConsistency() {
        VERIFY_MSG(segment_type_ != SegmentType::DiversitySegment and
                           chain_type_.Chain() != ImmuneChainType::AlphaTcrChain,
                   "TRA chain does not contain D gene segment");
        VERIFY_MSG(segment_type_ != SegmentType::DiversitySegment and
                           chain_type_.Chain() != ImmuneChainType::GammaTcrChain,
                   "TRG chain does not contain D gene segment");
        VERIFY_MSG(segment_type_ != SegmentType::DiversitySegment and
                           chain_type_.Chain() != ImmuneChainType::KappaIgChain,
                   "IGK chain does not contain D gene segment");
        VERIFY_MSG(segment_type_ != SegmentType::DiversitySegment and
                           chain_type_.Chain() != ImmuneChainType::LambdaIgChain,
                   "IGL chain does not contain D gene segment");
        return true;
    }

    std::ostream &operator<<(std::ostream &out, const ImmuneGeneType &gene_type) {
        out << gene_type.Chain() << gene_type.Segment();
        return out;
    }

    size_t ImmuneGeneTypeHasher::operator()(const ImmuneGeneType &gene_type) const {
        return std::hash<int>()(gene_type.Segment()) * ChainTypeHasher()(gene_type.Chain());
    }

}
