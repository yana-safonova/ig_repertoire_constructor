#pragma once

#include "verify.hpp"
#include "logger/logger.hpp"
#include "adt/concurrent_dsu.hpp"
#include "path_helper.hpp"

#include "sub_kmer_data.hpp"
#include "sub_kmer_splitter.hpp"
#include "sub_kmer_cutter.hpp"
#include "sub_kmer_block_file_reader.hpp"
#include "sub_kmer_block_file_writer.hpp"
#include "hamming_distance.hpp"

namespace ham_clustering {

template <class KMerData>
class ClusterizationIteration {
public:
    ClusterizationIteration() : next_iteration_(nullptr), block_threshold_(1e9), stride_(-1) {
    }

    ClusterizationIteration(std::shared_ptr <ClusterizationIteration> next_iteration, unsigned block_threshold, unsigned stride)
        : next_iteration_(next_iteration), block_threshold_(block_threshold), stride_(stride) {
    }

    void Run(const std::string & sub_kmers_file_name,
             const KMerData & data,
             ConcurrentDSU & components_dsu,
             bool need_logging) const;

private:
    void ProcessBlockQuadratic(ConcurrentDSU &components_dsu,
                               const std::vector<size_t> &block,
                               const KMerData &data) const;

    std::pair<size_t, size_t> SplitOnBlocks(const std::string & sub_kmers_file_name, const std::string & blocks_file_name) const;
    std::pair<size_t, size_t> ProcessAllBlocks(const KMerData & data,
                                               ConcurrentDSU & components_dsu,
                                               const std::string & blocks_file_name,
                                               const std::string & next_sub_kmers_file_name) const;

    std::shared_ptr <ClusterizationIteration> next_iteration_;
    unsigned block_threshold_;
    unsigned stride_;

    DECL_LOGGER("HammingClustering");
};

#if 1
static inline bool CanMerge(const ConcurrentDSU &, unsigned , unsigned ) {
    return true;
}
#else
#if 1
static bool CanMerge(const ConcurrentDSU &uf, unsigned x, unsigned y) {
    size_t szx = uf.set_size(x), szy = uf.set_size(y);
    const size_t hardthr = 2500;

    // Global threshold - no cluster larger than hard threshold
    if (szx + szy > hardthr)
        return false;

    // If one of the clusters is moderately large, than attach "almost" singletons
    // only.
    if ((szx > hardthr * 3 / 4 && szy > 50) ||
            (szy > hardthr * 3 / 4 && szx > 50))
        return false;

    return true;
}
#else
static bool CanMerge(const ConcurrentDSU &uf, unsigned x, unsigned y) {
    return (uf.set_size(x) + uf.set_size(y)) < 10000;
}
#endif
#endif

template <class KMerData>
void ClusterizationIteration<KMerData>::ProcessBlockQuadratic(ConcurrentDSU &components_dsu,
                                                       const std::vector<size_t> &block,
                                                       const KMerData &data) const {
    static unsigned tau = data.GetMaxCountMismatches();
    size_t blockSize = block.size();
    for (size_t i = 0; i < blockSize; ++i) {
        auto x = static_cast<unsigned>(block[i]);
        typename KMerData::KMer kmerx = data[x];

        for (size_t j = i + 1; j < blockSize; j++) {
            auto y = static_cast<unsigned>(block[j]);
            typename KMerData::KMer kmery = data[y];
            if (components_dsu.find_set(x) != components_dsu.find_set(y) &&
                CanMerge(components_dsu, x, y) &&
                GetHammingDistance <KMerData>(kmerx, kmery, tau) <= tau) {
                components_dsu.unite(x, y);
            }
        }
    }
}

template <class KMerData>
std::pair<size_t, size_t> ClusterizationIteration<KMerData>::SplitOnBlocks(const std::string & sub_kmers_file_name, const std::string & blocks_file_name) const {
    SubKMerSplitter <KMerData> splitter(sub_kmers_file_name, blocks_file_name);
    return splitter.split();
}

template <class KMerData>
std::pair<size_t, size_t> ClusterizationIteration<KMerData>::ProcessAllBlocks(const KMerData & data,
                                                           ConcurrentDSU & components_dsu,
                                                           const std::string & blocks_file_name,
                                                           const std::string & next_sub_kmers_file_name) const {
    SubKMerBlockFileReader <KMerData> blocks_reader(blocks_file_name, /* unlink */ true);

    SubKMerBlockFileWriter <KMerData> next_kmer_file(next_sub_kmers_file_name, data);
    size_t processed_blocks = 0;
    std::vector <size_t> block;
    size_t very_big_blocks_count = 0;
    while (blocks_reader.GetBlock(block)) {
        SubKMerCutterOnlyWithMismatchesFactory <KMerData> cutter_factory(data, block);
        if (block.size() > block_threshold_ && !cutter_factory.HasEmptyCutters()) {
            // Dump for next iteration.
            for (unsigned i = 0; i < cutter_factory.GetCuttersCount(); ++i) {
                next_kmer_file.Write(cutter_factory.Create(i), &block);
            }
            if (next_kmer_file.GetCountWrittenBlocks() % 10000 == 0) {
            	TRACE("Processed " << processed_blocks << " blocks, " << very_big_blocks_count << " from them are very big. "
                      "Produced " << next_kmer_file.GetCountWrittenBlocks() << " blocks for next iteration.");
            }
        } else {
            // Merge small blocks.
            if (block.size() >= 5000) {
                ++very_big_blocks_count;
            }
            ProcessBlockQuadratic(components_dsu, block, data);
        }
        ++processed_blocks;
        if (processed_blocks % 100000 == 0) {
            TRACE("Processed " << processed_blocks << " blocks, " << very_big_blocks_count << " from them are very big. "
	          "Produced " << next_kmer_file.GetCountWrittenBlocks() << " blocks for next iteration.");
        }
    }
    next_kmer_file.Close();

    return std::make_pair(processed_blocks, next_kmer_file.GetCountWrittenBlocks());
}

template <class KMerData>
void ClusterizationIteration<KMerData>::Run(const std::string & sub_kmers_file_name,
                                            const KMerData & data,
                                            ConcurrentDSU & components_dsu,
                                            bool need_logging) const {
    const std::string & blocks_file_name = sub_kmers_file_name + ".blocks";

    std::pair<size_t, size_t> stat = SplitOnBlocks(sub_kmers_file_name, blocks_file_name);
    if (need_logging) {
        TRACE("Splitting done."
             " Processed " << stat.first << " blocks."
             " Produced " << stat.second << " blocks.");
    }

    std::string next_sub_kmers_file_name = sub_kmers_file_name;
    stat = ProcessAllBlocks(data, components_dsu, blocks_file_name, next_sub_kmers_file_name);

    if (need_logging) {
        TRACE("Merge done."
             " Processed " << stat.first << " blocks."
             " Produced " << stat.second << " blocks for next iteration.");
    }
    if (stat.second != 0) {
        next_iteration_->Run(next_sub_kmers_file_name, data, components_dsu, need_logging);
    }
}

}
