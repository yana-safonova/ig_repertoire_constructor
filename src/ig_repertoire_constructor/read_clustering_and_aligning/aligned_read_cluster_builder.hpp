#pragma once

#include "aligned_read_cluster.hpp"
#include "read_archive.hpp"

#include "graph_walker.hpp"
#include "hamming_graph_building/kmer_data.hpp"

namespace ig_repertoire_constructor {

class AlignedReadClusterBuilder {
public:
    VectorAlignedReadClusterPtr Build(const ReadArchive & read_archive, const std::vector <std::vector <size_t> > & components) const;
};

VectorAlignedReadClusterPtr AlignedReadClusterBuilder::Build(const ReadArchive & read_archive, const std::vector <std::vector <size_t> > & components) const {
    KMerData kmer_data(read_archive, ig_cfg::get().aligning_params.min_overlap_length);
    INFO("KMers count: " << kmer_data.size());

    GraphWalker graph_walker(read_archive, kmer_data, components);
    VectorAlignedReadClusterPtr clusters_ptr(new std::vector <AlignedReadCluster>());
    for (size_t main_read_number = 0; main_read_number < read_archive.size(); ++main_read_number) {
        if (graph_walker.ReadIsUsed(main_read_number)) {
            continue;
        }
        AlignClusterBuilder clusterBuilder;
        graph_walker.Walk(main_read_number, clusterBuilder);
        clusters_ptr->push_back(clusterBuilder.GetCluster());
    }
    return clusters_ptr;
}

}
