#pragma once

#include "chain_type.hpp"

namespace germline_utils {

    enum SegmentType {
        UnknownSegment, VariableSegment, DiversitySegment, JoinSegment
    };

    std::ostream &operator<<(std::ostream &out, const SegmentType &segment_type);

    class ImmuneGeneType {
        ChainType chain_type_;
        SegmentType segment_type_;

        bool CheckChainSegmentConsistency();

    public:
        ImmuneGeneType() :
                chain_type_(),
                segment_type_(SegmentType::UnknownSegment) { }

        ImmuneGeneType(ChainType chain_type, SegmentType segment_type) :
                chain_type_(chain_type),
                segment_type_(segment_type) {
            CheckChainSegmentConsistency();
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