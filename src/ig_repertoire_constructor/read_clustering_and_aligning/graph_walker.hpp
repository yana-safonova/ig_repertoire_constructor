#pragma once

#include "verify.hpp"
#include "logger/logger.hpp"

#include "read_archive.hpp"

#include "hamming_graph_building/kmer_data.hpp"

namespace ig_repertoire_constructor {

class GraphWalker {
public:
    GraphWalker(const ReadArchive& read_archive,
                const KMerData & kmer_data,
                const std::vector <std::vector <size_t> > & components);

    bool ReadIsUsed(size_t read_number) const;
    void Walk(size_t main_read_number, AlignClusterBuilder & clusterBuilder);

private:
    void VisitAllFromKMer(std::vector <size_t> & queue, size_t read_number, size_t kmer_number);
    void TryVisitAnotherKmer(std::vector <size_t> & queue, const KMerInfo & kmer_info, const KMerInfo & another_kmer_info);

    const KMerData & kmer_data_;
    const std::vector <std::vector <size_t> > & components_;
    std::vector <bool> read_is_used_;
    std::vector <bool> component_is_used_;
    std::vector <int> shift_for_read_;
    std::vector <size_t> component_number_by_kmer_;
    const unsigned threshold_for_single_shift_;

    DECL_LOGGER("GraphWalker");
};

inline GraphWalker::GraphWalker(const ReadArchive& read_archive,
            const KMerData & kmer_data,
            const std::vector <std::vector <size_t> > & components)
            : kmer_data_(kmer_data), components_(components),
              threshold_for_single_shift_(ig_cfg::get().aligning_params.threshold_for_single_shift) {

    read_is_used_.assign(read_archive.size(), false);
    component_is_used_.assign(components.size(), false);
    shift_for_read_.assign(read_archive.size(), 0);

    component_number_by_kmer_.assign(kmer_data.size(), -1);
    for (size_t component_number = 0; component_number < components_.size(); ++component_number) {
        for (const auto & vertex : components_[component_number]) {
            component_number_by_kmer_[vertex] = component_number;
        }
    }
}

inline bool GraphWalker::ReadIsUsed(size_t read_number) const {
    return read_is_used_[read_number];
}

inline void GraphWalker::Walk(size_t main_read_number, AlignClusterBuilder & clusterBuilder) {
    DEBUG("Start walking from read #" << main_read_number);
    VERIFY(!ReadIsUsed(main_read_number));

    read_is_used_[main_read_number] = true;
    shift_for_read_[main_read_number] = 0;
    std::vector <size_t> queue;
    queue.push_back(main_read_number);

    for (size_t head = 0; head < queue.size(); ++head) {
        size_t read_number = queue[head];
        clusterBuilder.Add(read_number, shift_for_read_[read_number]);

        std::pair <size_t, size_t> range = kmer_data_.GetRange(read_number);
        for (size_t kmer_number = range.first; kmer_number < range.second; ++kmer_number) {
            VisitAllFromKMer(queue, read_number, kmer_number);
        }
    }
}

inline void GraphWalker::VisitAllFromKMer(std::vector <size_t> & queue, size_t read_number, size_t kmer_number) {
    size_t component_number = component_number_by_kmer_[kmer_number];
    if (component_is_used_[component_number]) {
        return;
    }

    DEBUG("Visit component @" << component_number << " from kmer $" << kmer_number << " of read #" << read_number);

    component_is_used_[component_number] = true;
    const KMerInfo & kmer_info = kmer_data_.GetInfo(kmer_number);
    VERIFY(kmer_info.read_number == read_number);
    for (size_t another_kmer_number : components_[component_number]) {
        const KMerInfo & another_kmer_info = kmer_data_.GetInfo(another_kmer_number);
        TryVisitAnotherKmer(queue, kmer_info, another_kmer_info);
    }
}

inline void GraphWalker::TryVisitAnotherKmer(std::vector <size_t> & queue, const KMerInfo & kmer_info, const KMerInfo & another_kmer_info) {
    size_t read_number = kmer_info.read_number;
    size_t another_read_number = another_kmer_info.read_number;
    int new_shift_for_read = shift_for_read_[read_number] + kmer_info.offset - another_kmer_info.offset;
    if (abs(size_t(new_shift_for_read)) >= threshold_for_single_shift_) {
    	return;
    }
    if (read_is_used_[another_read_number]) {
        if (shift_for_read_[another_read_number] != new_shift_for_read) {
            WARN("For read #" << another_read_number << " found two optimal shifts: " << shift_for_read_[another_read_number] << " and " << new_shift_for_read << " from read #" << read_number);
        }
        return;
    }
    DEBUG("Found new shift(" << new_shift_for_read << ") for read #" << another_read_number << " using read #" << read_number);
    read_is_used_[another_read_number] = true;
    shift_for_read_[another_read_number] = new_shift_for_read;
    queue.push_back(another_read_number);
}

}
