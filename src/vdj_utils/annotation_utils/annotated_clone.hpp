#pragma once

#include "cdr_labeling_primitives.hpp"
#include <read_archive.hpp>
#include <annotation_utils/aa_annotation/aa_annotation.hpp>
#include <annotation_utils/shm_annotation/shm_annotation.hpp>
#include "../alignment_utils/pairwise_alignment.hpp"

namespace annotation_utils {
    enum StructuralRegion { CDR, FR, CDR1, CDR2, CDR3, FR1, FR2, FR3, FR4, UnknownRegion, AnyRegion };

    std::ostream& operator<<(std::ostream& out, const StructuralRegion &region);

    class AnnotatedClone {
        core::Read read_;

        std::map<StructuralRegion, seqan::Dna5String> region_string_map_;
        std::map<StructuralRegion, CDRRange> region_range_map_;

        alignment_utils::ImmuneGeneReadAlignment v_alignment_;
        alignment_utils::ImmuneGeneReadAlignment j_alignment_;

        AminoAcidAnnotation<core::Read> aa_annotation_;

        GeneSegmentSHMs v_shms_;
        GeneSegmentSHMs j_shms_;

        size_t size_;

        void CheckRangeConsistencyFatal(CDRRange range);

        void UpdateStructuralRegion(StructuralRegion region, CDRRange range);

        void Initialize(CDRLabeling cdr_labeling);

    public:
        AnnotatedClone(core::Read read,
                       CDRLabeling cdr_labeling,
                       alignment_utils::ImmuneGeneReadAlignment v_alignment,
                       alignment_utils::ImmuneGeneReadAlignment j_alignment,
                       AminoAcidAnnotation<core::Read> aa_annotation,
                       GeneSegmentSHMs v_shms,
                       GeneSegmentSHMs j_shms) :
                read_(read),
                v_alignment_(v_alignment),
                j_alignment_(j_alignment),
                aa_annotation_(aa_annotation),
                v_shms_(v_shms),
                j_shms_(j_shms),
                size_(1) {
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


        CDRRange CDR1Range() const {
            return GetRangeByRegion(StructuralRegion::CDR1);
        }

        CDRRange CDR2Range() const {
            return GetRangeByRegion(StructuralRegion::CDR2);
        }

        CDRRange CDR3Range() const {
            return GetRangeByRegion(StructuralRegion::CDR3);
        }

        const core::Read& Read() const { return read_; }

        const alignment_utils::ImmuneGeneReadAlignment& VAlignment() const { return v_alignment_; }

        const alignment_utils::ImmuneGeneReadAlignment& JAlignment() const { return j_alignment_; }

        germline_utils::ChainType ChainType() const { return v_alignment_.subject().Chain(); }

        const GeneSegmentSHMs& VSHMs() const { return v_shms_; }

        const GeneSegmentSHMs& JSHMs() const { return j_shms_; }

        bool HasStopCodon() const { return aa_annotation_.HasStopCodon(); }

        bool InFrame() const { return aa_annotation_.InFrame(); }

        bool Productive() const { return !HasStopCodon() and InFrame(); }

        seqan::String<seqan::AminoAcid> AA() const { return aa_annotation_.AA(); }

        const germline_utils::ImmuneGene& VGene() const { return v_alignment_.subject(); }

        const germline_utils::ImmuneGene& JGene() const { return j_alignment_.subject(); }

        size_t ORF() const { return aa_annotation_.ORF(); }

        size_t  Size() const { return size_; }

        void SetSize(size_t n) { size_ = n; }

        const alignment_utils::ImmuneGeneReadAlignment& GetAlignmentBySegment(
                germline_utils::SegmentType segment_type) const;

        size_t GetAminoAcidPosByNucleotidePos(size_t nucl_pos) const;

        char GetAminoAcidByNucleotidePos(size_t nucl_pos) const;

        StructuralRegion GetRegionBySHM(SHM shm) const;

        seqan::Dna5String GetJNucleotides() const;

        seqan::Dna5String GetCDR3JNucleotides() const;

    };

    std::ostream& operator<<(std::ostream& out, const AnnotatedClone &obj);
}