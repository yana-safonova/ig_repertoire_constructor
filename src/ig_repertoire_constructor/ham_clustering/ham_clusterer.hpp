#pragma once

#include "verify.hpp"
#include "logger/logger.hpp"
#include "adt/concurrent_dsu.hpp"
#include "path_helper.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "sub_kmer_cutter.hpp"
#include "sub_kmer_block_file_writer.hpp"
#include "clusterization_iteration.hpp"
#include "hamming_distance.hpp"

namespace ham_clustering {

template <class KMerData, class SubKMerCutter>
class KMerHamClustererImpl {
public:
    KMerHamClustererImpl(const std::vector <SubKMerCutter> & sub_kmer_cutters, std::shared_ptr <ClusterizationIteration <KMerData> > first_iteration)
        : sub_kmer_cutters_(sub_kmer_cutters), first_iteration_(first_iteration) {
    }

    void Run(const KMerData & data, ConcurrentDSU & components_dsu, const std::string & prefix, bool need_logging) const {
        const std::string & first_sub_kmers_file_name = prefix + "ham_clustering.sub_kmers";
        SubKMerBlockFileWriter <KMerData> first_kmer_file(first_sub_kmers_file_name, data);
        for (const auto & cutter : sub_kmer_cutters_) {
            first_kmer_file.Write(cutter, NULL);
        }
        first_kmer_file.Close();
        first_iteration_->Run(first_sub_kmers_file_name, data, components_dsu, need_logging);
    }

private:
    const std::vector <SubKMerCutter> & sub_kmer_cutters_;
    std::shared_ptr <ClusterizationIteration <KMerData> > first_iteration_;
};

template <class KMerData>
class KMerHamClusterer {
public:
    KMerHamClusterer() {}

    void Cluster(const std::string &prefix, const KMerData &data, std::vector < std::vector <size_t> > & components, unsigned threads_count) const;
private:
    inline void RunImpl(unsigned tau, const KMerData &data, ConcurrentDSU & components_dsu, const std::string & prefix,
                        const std::vector <SubKMerCutterBySubPermutation <KMerData> > & cutters) const;

    DECL_LOGGER("HammingClustering");
};

template <class KMerData>
inline void KMerHamClusterer<KMerData>::RunImpl(
                    unsigned , const KMerData &data, ConcurrentDSU & components_dsu, const std::string & prefix,
                    const std::vector <SubKMerCutterBySubPermutation <KMerData> > & cutters) const {
    bool need_logging = (data.size() > 10000);
    auto third_iteration = std::make_shared <ClusterizationIteration <KMerData> >();
    auto second_iteration = std::make_shared <ClusterizationIteration <KMerData> >(third_iteration, 15000, 0);
    auto first_iteration = std::make_shared <ClusterizationIteration <KMerData> >(second_iteration, 50, 0);
    KMerHamClustererImpl <KMerData, SubKMerCutterBySubPermutation <KMerData> > runnable(cutters, first_iteration);
    runnable.Run(data, components_dsu, prefix, need_logging);
}

template <class KMerData>
inline void KMerHamClusterer<KMerData>::Cluster(const std::string & prefix, const KMerData & data, std::vector < std::vector <size_t> > & components, unsigned threads_count) const {

    ConcurrentDSU components_dsu(data.size());

    unsigned tau = data.GetMaxCountMismatches();
    if (data.size() <= 100) {
        for (unsigned i = 0; i < data.size(); ++i) {
            for (unsigned j = i + 1; j < data.size(); ++j) {
                if (components_dsu.find_set(i) != components_dsu.find_set(j) &&
                    GetHammingDistance <KMerData>(data[i], data[j], tau) <= tau) {
                    components_dsu.unite(i, j);
                }
            }
        }
    } else {
        SubKMerCutterWithStrideFactory <KMerData> cutter_factory(data.GetKmerLength(), tau + 1, tau);

    #pragma omp parallel for num_threads(std::min(tau + 1, threads_count))
        for (unsigned i = 0; i <= tau; ++i) {
            auto cutter = cutter_factory.Create(i);
            std::vector <SubKMerCutterBySubPermutation <KMerData> > cutters(1, cutter);
            RunImpl(tau, data, components_dsu, prefix + char('a' + i) + "_", cutters);
        }
    }

    components_dsu.get_sets(components);
}

}
