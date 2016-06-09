#pragma once

#include "cdr_labeling_primitives.hpp"
#include <read_archive.hpp>
#include "../alignment_utils/pairwise_alignment.hpp"

namespace annotation_utils {
    enum StructuralRegion { CDR1, CDR2, CDR3 };

    std::ostream& operator<<(std::ostream& out, const StructuralRegion &region);

    class CDRAnnotatedClone {
        const core::Read &read_;

        std::unordered_map<StructuralRegion, seqan::Dna5String, std::hash<int>> region_string_map_;
        std::unordered_map<StructuralRegion, CDRRange, std::hash<int>> region_range_map_;

        alignment_utils::ImmuneGeneReadAlignment v_alignment_;
        alignment_utils::ImmuneGeneReadAlignment j_alignment_;

        void CheckRangeConsistencyFatal(CDRRange range);

        void UpdateStructuralRegion(StructuralRegion region, CDRRange range);

        void Initialize(CDRLabeling cdr_labeling);

    public:
        CDRAnnotatedClone(const core::Read &read,
                          CDRLabeling cdr_labeling,
                          alignment_utils::ImmuneGeneReadAlignment v_alignment,
                          alignment_utils::ImmuneGeneReadAlignment j_alignment) : read_(read),
                                                                                  v_alignment_(v_alignment),
                                                                                  j_alignment_(j_alignment) {
            Initialize(cdr_labeling);
        }

        bool RegionIsEmpty(StructuralRegion region) const;

        CDRRange GetRangeByRegion(StructuralRegion region) const;

        seqan::Dna5String GetRegionString(StructuralRegion region) const;

        seqan::Dna5String CDR1() const {
            return GetRegionString(StructuralRegion::CDR1);
        }

        seqan::Dna5String CDR2() const {
            return GetRegionString(StructuralRegion::CDR2);
        }

        seqan::Dna5String CDR3() const {
            return GetRegionString(StructuralRegion::CDR3);
        }

        const core::Read& Read() const { return read_; }

        const alignment_utils::ImmuneGeneReadAlignment& VAlignment() const { return v_alignment_; }

        const alignment_utils::ImmuneGeneReadAlignment& JAlignment() const { return j_alignment_; }
    };

    std::ostream& operator<<(std::ostream& out, const CDRAnnotatedClone &obj);
}