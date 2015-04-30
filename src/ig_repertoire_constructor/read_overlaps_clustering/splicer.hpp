#pragma once

#include <memory>
#include "splicer_result.hpp"
#include "overlap_borders_finder.hpp"

namespace ig_repertoire_constructor {

class Splicer {
public:
    std::shared_ptr <SplicerResult> SpliceOverlap(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                  const AlignedReadCluster & cluster,
                                                  const std::vector <size_t> & indices) const;
};

inline std::shared_ptr <SplicerResult> Splicer::SpliceOverlap(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                              const AlignedReadCluster & cluster,
                                                              const std::vector <size_t> & indices) const {
    auto bordersFinder = OverlapBordersFinderFactory().Create();
    auto overlap_borders = bordersFinder->FindOptimalOverlapBorders(read_archive_ptr, cluster, indices);

    auto result = std::make_shared <SplicerResult>();
    result->spliced_reads.reserve(indices.size());
    for (auto index : indices) {
        const auto & aligned_read = cluster[index];
        const io::SingleRead & read = read_archive_ptr->operator[](aligned_read.read_number);
        if (aligned_read.shift <= overlap_borders.first
                && overlap_borders.first <= overlap_borders.second
                && overlap_borders.second <= (unsigned)read.size() + aligned_read.shift) {
            unsigned from = overlap_borders.first - aligned_read.shift;
            unsigned to = overlap_borders.second - aligned_read.shift;
            result->spliced_reads.emplace_back(read_archive_ptr, aligned_read.read_number, from, to);
            result->spliced_read_indices.push_back(index);
        } else {
            result->bad_read_numbers.push_back(aligned_read.read_number);
        }
    }
    return result;
}

}
