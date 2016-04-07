#pragma once

#include "chain_type.hpp"

#include <vector>

namespace germline_utils {
    enum SegmentType {
        UnknownSegment, VariableSegment, DiversitySegment, JoinSegment
    };

    std::ostream &operator<<(std::ostream &out, const SegmentType &segment_type);

    struct SegmentTypeHasher {
        size_t operator()(const SegmentType &segment_type) const { return std::hash<int>()(int(segment_type)); }
    };

    class ImmuneGeneType {
        ChainType chain_type_;
        SegmentType segment_type_;

        bool CheckChainSegmentConsistency();

        void Initialize(ChainType chain_type, SegmentType segment_type);

    public:
        ImmuneGeneType() :
                chain_type_(),
                segment_type_(SegmentType::UnknownSegment) { }

        ImmuneGeneType(std::string immune_gene_str);

        ImmuneGeneType(ChainType chain_type, SegmentType segment_type) {
            Initialize(chain_type, segment_type);
        }

        LymphocyteType Lymphocyte() const { return chain_type_.Lymphocyte(); }

        ChainType Chain() const { return chain_type_; }

        SegmentType Segment() const { return segment_type_; }

        bool operator==(const ImmuneGeneType &gene_type) const {
            return chain_type_ == gene_type.Chain() and segment_type_ == gene_type.Segment();
        }
    };

    std::ostream &operator<<(std::ostream &out, const ImmuneGeneType &gene_type);

    struct ImmuneGeneTypeHasher {
        size_t operator()(const ImmuneGeneType& gene_type) const;
    };
}