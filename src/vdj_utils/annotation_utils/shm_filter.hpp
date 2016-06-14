#pragma once

#include "shm_annotation.hpp"

namespace annotation_utils {
    class SHMFilter {
    public:
        virtual bool FilterSHM(SHM shm) const = 0; // return true if SHM should be fitered out

        virtual ~SHMFilter() { }
    };

    class PositionalSHMFilter {
        const GeneSegmentSHMs &shms_;
        size_t max_num_skipped_start_;
        size_t max_num_skipped_end_;

        size_t first_meaning_read_pos_; // first position on read corresponding to good SHMs
        size_t first_meaning_gene_pos_; // first position on gene corresponding to good SHMs
        size_t last_meaning_read_pos_; // last position on read corresponding to good SHMs
        size_t last_meaning_gene_pos_; // last position on gene corresponding to good SHMs

        void ComputeStartMeaningPositions();

        void ComputeEndMeaningPositions();

        void ComputeMeaningPositions();

    public:
        PositionalSHMFilter(const GeneSegmentSHMs &shms,
                            size_t max_num_skipped_start,
                            size_t max_num_skipped_end) : shms_(shms),
                                                          max_num_skipped_start_(max_num_skipped_start),
                                                          max_num_skipped_end_(max_num_skipped_end) {
            ComputeMeaningPositions();
        }

        bool FilterSHM(SHM shm) const;
    };
}