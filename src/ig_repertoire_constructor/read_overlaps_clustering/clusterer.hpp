#pragma once

#include "sequence/rtseq.hpp"
#include "ham_clustering/ham_clusterer.hpp"
#include "kmer_data.hpp"

namespace ig_repertoire_constructor {

class Clusterer {
public:
    void ProcessCluster(std::shared_ptr <ReadArchive> read_archive_ptr,
                        const AlignedReadCluster & aligned_read_cluster,
                        std::vector <std::vector <SplicedRead> > & components,
                        const std::string & prefix) const;

private:
    void MakeHammingClusterization(const std::vector <SplicedRead> & spliced_reads,
                                   std::vector <std::vector <size_t> > & clusters,
                                   const std::string & prefix) const;
};

inline void Clusterer::ProcessCluster(std::shared_ptr <ReadArchive> read_archive_ptr,
                                      const AlignedReadCluster & aligned_read_cluster,
                                      std::vector <std::vector <SplicedRead> > & components,
                                      const std::string & prefix) const {
    components.clear();

    std::vector <std::vector <size_t> > clusters;

    std::vector <size_t> all_indices(aligned_read_cluster.size());
    for (size_t index = 0; index < all_indices.size(); ++index) {
        all_indices[index] = index;
    }
    std::vector <std::vector <size_t> > stack;
    stack.push_back(all_indices);

    while (!stack.empty()) {
        std::vector <size_t> indices = stack.back();
        stack.pop_back();

        auto splicer_result = Splicer().SpliceOverlap(read_archive_ptr, aligned_read_cluster, indices);
        for (size_t bad_read_number : splicer_result->bad_read_numbers) {
            const auto & read = read_archive_ptr->operator [](bad_read_number);
            SplicedRead spliced_read(read_archive_ptr, bad_read_number, 0, (unsigned)read.size());
            components.push_back(std::vector <SplicedRead>(1, spliced_read));
        }

        MakeHammingClusterization(splicer_result->spliced_reads, clusters, prefix);
        if (clusters.size() == 1) {
            components.push_back(splicer_result->spliced_reads);
        } else {
            std::vector <size_t> sub_indices;
            for (const auto & cluster : clusters) {
                sub_indices.resize(cluster.size());
                for (size_t i = 0; i < cluster.size(); ++i) {
                    sub_indices[i] = splicer_result->spliced_read_indices[cluster[i]];
                }
                stack.push_back(sub_indices);
            }
        }
    }
}

inline void Clusterer::MakeHammingClusterization(const std::vector <SplicedRead> & spliced_reads,
                                                 std::vector <std::vector <size_t> > & clusters,
                                                 const std::string & prefix) const {
    SplicedKMerData kmer_data(spliced_reads, ig_cfg::get().aligning_params.overlap_mismatches_threshold);
    ham_clustering::KMerHamClusterer <SplicedKMerData> ham_clusterer;
    ham_clusterer.Cluster(prefix, kmer_data, clusters, 1); // NOTE There is in multi-thread mode yet.
}

}
