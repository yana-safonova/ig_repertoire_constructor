#pragma once

#include <seqan/seq_io.h>
#include "utils.hpp"
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"
#include "umi_utils.hpp"

namespace clusterer {

    typedef std::function<size_t(const seqan::Dna5String &, const seqan::Dna5String &)> DistFunction;

    struct ClusteringMode {
        const DistFunction dist;
        const size_t threshold;

        ClusteringMode(DistFunction dist_, size_t threshold_) : dist(dist_), threshold(threshold_) { }

        static const ClusteringMode hamming;
        static const ClusteringMode edit;
    };


    class ReflexiveUmiPairsIterator {
    public:
        ReflexiveUmiPairsIterator(size_t current, size_t last) : current_(current), last_(last) {}

        ReflexiveUmiPairsIterator operator ++();
        ReflexiveUmiPairsIterator operator ++(int);
        bool operator==(ReflexiveUmiPairsIterator other) const;
        bool operator!=(ReflexiveUmiPairsIterator other) const;
        std::pair<size_t, size_t> operator*() const;

    private:
        size_t current_;
        size_t last_;
    };

    // TODO: add GraphUmiPairsIterable/Iterator and probably full graph

    class ReflexiveUmiPairsIterable {
    public:
        ReflexiveUmiPairsIterable(size_t count): count_(count) {}

        ReflexiveUmiPairsIterator begin() const;
        ReflexiveUmiPairsIterator end() const;

    private:
        size_t count_;
    };


    struct SingleReadClusterElement {
        SingleReadClusterElement(const seqan::Dna5String& read) : center(read), weight(1) {}

        seqan::Dna5String center;
        size_t weight;
    };

    template <typename ElementType>
    struct Cluster {
    public:
        Cluster(ElementType element_ = ElementType(), size_t weight_ = 1) : center(element_), weight(weight_), members{element_} {}

        seqan::Dna5String center;
        size_t weight;
        std::unordered_set<ElementType> members;
    };

    template <typename ElementType>
    using ClusterPtr = std::shared_ptr<Cluster<ElementType>>;


    template <typename ElementType, typename UmiPairsIterable>
    class Clusterer {
    public:
        static ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>> cluster(
                const ClusteringMode& mode,
                const std::vector<UmiPtr>& umis,
                const ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>>& umis_to_clusters,
                const UmiPairsIterable& umi_pairs_iterable);

    private:
        static ClusterPtr<ElementType> mergeClusters(const ClusterPtr<ElementType> first, const ClusterPtr<ElementType> second);
        static seqan::Dna5String findNewCenter(std::unordered_set<ElementType>& members);
    };

    template <typename ElementType, typename UmiPairsIterable>
    ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>> Clusterer<ElementType, UmiPairsIterable>::cluster(
            const ClusteringMode& mode,
            const std::vector<UmiPtr>& umis,
            const ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>>& umis_to_clusters,
            const UmiPairsIterable& umi_pairs_iterable) {
        ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>> result(umis_to_clusters);
        for (auto& umi_pair : umi_pairs_iterable) {
            UmiPtr first_umi = umis[umi_pair.first];
            UmiPtr second_umi = umis[umi_pair.second];
            for (auto& first_cluster : umis_to_clusters.forth(first_umi)) {
                for (auto& second_cluster : umis_to_clusters.forth(second_umi)) {
                    if (first_cluster == second_cluster) continue;
                    if (mode.dist(first_cluster->getRepresentative(), second_cluster->getRepresentative()) <= mode.threshold) {
                        VERIFY_MSG(result.removeSecond(first_cluster), "Trying to remove an absent cluster");
                        VERIFY_MSG(result.removeSecond(second_cluster), "Trying to remove an absent cluster");
                        auto merged_cluster = mergeClusters(first_cluster, second_cluster);
                        // TODO: avoid returning copied umi set by providing access to its begin() and end()
                        const auto& first_cluster_umis = umis_to_clusters.back(first_cluster);
                        const auto& second_cluster_umis = umis_to_clusters.back(second_cluster);
                        std::unordered_set<UmiPtr> merged_umis;
                        std::set_union(first_cluster_umis.begin(), first_cluster_umis.end(),
                                       second_cluster_umis.begin(), second_cluster_umis.end(),
                                       merged_umis);
                        result.add(merged_umis, merged_cluster);
                    }
                }
            }
        }
        return result;
    };

    template <typename ElementType, typename UmiPairsIterable>
    ClusterPtr<ElementType> Clusterer<ElementType, UmiPairsIterable>::mergeClusters(const ClusterPtr<ElementType> first, const ClusterPtr<ElementType> second) {
        std::unordered_set<ElementType> members;
        std::set_union(first->members.begin(), first->members.end(),
                       second->members.begin(), second->members.end(),
                       members);
        VERIFY_MSG(members.size() == first->members.size() + second->members.size(), "Clusters already intersect.");
        if (first->members.size() < second->members.size()) {
            std::swap(first, second);
        }
        size_t max_size = first->members.size();
        bool need_new_center = max_size <= 10 || __builtin_clzll(max_size) != __builtin_clz(members.size());
        ElementType center = need_new_center ? findNewCenter(members) : first->center;
        return ClusterPtr<ElementType>(center, first->weight + second->weight, members);
    }

    template <typename ElementType, typename UmiPairsIterable>
    seqan::Dna5String Clusterer<ElementType, UmiPairsIterable>::findNewCenter(std::unordered_set<ElementType>& members) {
        VERIFY_MSG(members.size() >= 2, "Too few elements to find new center.");
        seqan::Dna5String center = members.begin()->center;
        std::vector<std::vector<size_t> > cnt(length(center), std::vector<size_t>(4));
        for (auto& member : members) {
            size_t limit = std::min(length(member->center), length(center));
            for (size_t i = 0; i < limit; i ++) {
                cnt[i][member->center[i]] ++;
            }
        }
        for (size_t i = 0; i < length(center); i ++) {
            for (size_t j = 0; j < cnt[i].size(); j ++) {
                if (cnt[i][j] > cnt[i][center[i]]) {
                    center[i] = j;
                }
            }
        }
        return center;
    }
}