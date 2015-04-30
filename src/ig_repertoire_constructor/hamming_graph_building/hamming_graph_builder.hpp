#pragma once

#include "aligned_read_cluster.hpp"
#include "read_archive.hpp"

#include "kmer_data.hpp"
#include "ham_clustering/ham_clusterer.hpp"

namespace ig_repertoire_constructor {

class HammingGraphBuilder {
public:
    std::shared_ptr <std::vector <std::vector <size_t> > > Build(const ReadArchive & read_archive) const;
};

std::shared_ptr <std::vector <std::vector <size_t> > > HammingGraphBuilder::Build(const ReadArchive & read_archive) const {
    KMerData kmer_data(read_archive, ig_cfg::get().aligning_params.min_overlap_length);
    INFO("KMers count: " << kmer_data.size());

    ham_clustering::KMerHamClusterer <KMerData> clusterer;
    auto components_ptr = make_shared <std::vector <std::vector <size_t> > >();
    std::string prefix = path::append_path(ig_cfg::get().io.temp_files, "temp_");
    clusterer.Cluster(prefix, kmer_data, *components_ptr, ig_cfg::get().rp.threads_count);
    return components_ptr;
}

}
