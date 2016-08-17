#pragma once

#include <alignment_utils/pairwise_alignment.hpp>
#include <annotation_utils/cdr_labeling_primitives.hpp>
#include "shm_annotation.hpp"
#include "../aa_annotation/aa_annotation.hpp"

namespace annotation_utils {
    class BaseSHMCalculator {
    public:
        virtual GeneSegmentSHMs ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                            const AminoAcidAnnotation<core::Read>& aa_annotation,
                                            const CDRLabeling &cdr_labeling) = 0;

        virtual ~BaseSHMCalculator() { }
    };

    // NaiveSHMCalculator interprets all differences in alignment between Gene and Read as somatic hypermutations
    class NaiveSHMCalculator : public BaseSHMCalculator {
    public:
        GeneSegmentSHMs ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                    const AminoAcidAnnotation<core::Read>& aa_annotation,
                                    const CDRLabeling &cdr_labeling);
    };

    // StartEndFilteringSHMCalculator ignores differences at the start and end of alignment since
    // they typically present uncorrected errors and gene cleavage
    class StartEndFilteringSHMCalculator : public BaseSHMCalculator {
        size_t max_skipped_start_;
        size_t max_skipped_end_;

        size_t first_meaning_read_pos_; // first position on read corresponding to good SHMs
        size_t first_meaning_gene_pos_; // first position on gene corresponding to good SHMs
        size_t last_meaning_read_pos_; // last position on read corresponding to good SHMs
        size_t last_meaning_gene_pos_; // last position on gene corresponding to good SHMs

        void ComputeStartMeaningPositions(const GeneSegmentSHMs &all_shms);

        void ComputeEndMeaningPositions(const GeneSegmentSHMs &all_shms,
                                        size_t end_read_pos, size_t end_gene_pos);

        void ComputeMeaningPositions(const GeneSegmentSHMs& all_shms,
                                     const alignment_utils::ImmuneGeneReadAlignment& alignment);

    public:
        StartEndFilteringSHMCalculator(size_t max_skipped_start, size_t max_skipped_end) :
                max_skipped_start_(max_skipped_start),
                max_skipped_end_(max_skipped_end) { }

        GeneSegmentSHMs ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                    const AminoAcidAnnotation<core::Read>& aa_annotation,
                                    const CDRLabeling &cdr_labeling);
    };

    // CDRFilteringSHMCalculator ignores SHMs in CDR3
    class CDRFilteringSHMCalculator : public BaseSHMCalculator {
        size_t first_meaning_read_pos_; // first position on read corresponding to good SHMs
        size_t first_meaning_gene_pos_; // first position on gene corresponding to good SHMs
        size_t last_meaning_read_pos_; // last position on read corresponding to good SHMs
        size_t last_meaning_gene_pos_; // last position on gene corresponding to good SHMs

        void ComputeStartMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                          const GeneSegmentSHMs& all_shms,
                                          const CDRLabeling &cdr_labeling);

        void ComputeEndMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                        const GeneSegmentSHMs& all_shms,
                                        const CDRLabeling &cdr_labeling);

        void ComputeMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                     const GeneSegmentSHMs& all_shms,
                                     const CDRLabeling &cdr_labeling);

    public:
        GeneSegmentSHMs ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                    const AminoAcidAnnotation<core::Read>& aa_annotation,
                                    const CDRLabeling &cdr_labeling);
    };
}