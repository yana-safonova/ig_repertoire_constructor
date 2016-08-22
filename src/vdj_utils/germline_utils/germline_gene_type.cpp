#include <verify.hpp>

#include "germline_gene_type.hpp"

#include <vector>

namespace germline_utils {

    std::vector<std::string> segment_type_strs { "Unknown gene segment", "V", "D", "J" };

    std::ostream &operator<<(std::ostream &out, const SegmentType &segment_type) {
        out << segment_type_strs[int(segment_type)];
        return out;
    }

    void check_segment_str_correctness_fatal(std::string segment_str) {
        bool segment_str_correct = false;
        for(auto it = segment_type_strs.begin(); it != segment_type_strs.end(); it++)
            if(*it == segment_str)
                return;
        VERIFY_MSG(segment_str_correct, "Segment type is unknown: " << segment_str);
    }

    SegmentType convert_set_to_segment_type(std::string segment_str) {
        for(size_t i = 0; i < segment_type_strs.size(); i++)
            if(segment_type_strs[i] == segment_str)
                return SegmentType(i);
        return SegmentType::UnknownSegment;
    }

    void ImmuneGeneType::Initialize(ChainType chain_type, SegmentType segment_type) {
        chain_type_ = chain_type;
        segment_type_ = segment_type;
        CheckChainSegmentConsistency();
    }

    ImmuneGeneType::ImmuneGeneType(std::string immune_gene_str) {
        VERIFY_MSG(immune_gene_str.size() == 4, "Format of immune gene string is not correct: " << immune_gene_str);
        std::string segment_str = immune_gene_str.substr(3, 1);
        check_segment_str_correctness_fatal(segment_str);
        SegmentType segment_type = convert_set_to_segment_type(segment_str);
        Initialize(ChainType(immune_gene_str.substr(0, 3)), segment_type);
    }

    bool ImmuneGeneType::CheckChainSegmentConsistency() {
        VERIFY_MSG(!(segment_type_ == SegmentType::DiversitySegment and
                     chain_type_.Chain() == ImmuneChainType::AlphaTcrChain),
                   "TRA chain does not contain D gene segment");
        VERIFY_MSG(!(segment_type_ == SegmentType::DiversitySegment and
                     chain_type_.Chain() == ImmuneChainType::GammaTcrChain),
                   "TRG chain does not contain D gene segment");
        VERIFY_MSG(!(segment_type_ == SegmentType::DiversitySegment and
                     chain_type_.Chain() == ImmuneChainType::KappaIgChain),
                   "IGK chain does not contain D gene segment");
        VERIFY_MSG(!(segment_type_ == SegmentType::DiversitySegment and
                     chain_type_.Chain() == ImmuneChainType::LambdaIgChain),
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
