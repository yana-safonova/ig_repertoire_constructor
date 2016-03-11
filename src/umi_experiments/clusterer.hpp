#pragma once

#include <seqan/seq_io.h>
#include "utils.hpp"
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"
#include "../graph_utils/sparse_graph.hpp"
#include "umi_utils.hpp"
#include "disjoint_sets.hpp"

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

    class ReflexiveUmiPairsIterable {
    public:
        ReflexiveUmiPairsIterable(size_t count): count_(count) {}

        ReflexiveUmiPairsIterator begin() const;
        ReflexiveUmiPairsIterator end() const;

    private:
        size_t count_;
    };

    class GraphUmiPairsIterator {
    public:
        GraphUmiPairsIterator(const GraphUmiPairsIterator& other) = default;
        GraphUmiPairsIterator(const SparseGraphPtr& graph, size_t vertex, const SparseGraph::EdgesIterator& current_edge)
                : graph_(graph), vertex_(vertex), current_edge_(std::make_shared<SparseGraph::EdgesIterator>(current_edge)), advances_(0) {}

        GraphUmiPairsIterator operator ++();
        GraphUmiPairsIterator operator ++(int);
        bool operator==(GraphUmiPairsIterator other) const;
        bool operator!=(GraphUmiPairsIterator other) const;
        std::pair<size_t, size_t> operator*() const;

    private:
        const SparseGraphPtr graph_;
        size_t vertex_;
        std::shared_ptr<SparseGraph::EdgesIterator> current_edge_;
        // just for asserts
        size_t advances_;

        void advance();
    };

    class GraphUmiPairsIterable {
    public:
        GraphUmiPairsIterable(const SparseGraphPtr& graph);

        GraphUmiPairsIterator begin() const;
        GraphUmiPairsIterator end() const;

    private:
        const SparseGraphPtr& graph_;
        size_t first_connected_vertex_;
    };

    class FullGraphUmiPairsIterator {
    public:
        FullGraphUmiPairsIterator(size_t current_first, size_t current_second, size_t last)
                : current_first_(current_first), current_second_(current_second), last_(last) {}

        FullGraphUmiPairsIterator operator ++();
        FullGraphUmiPairsIterator operator ++(int);
        bool operator==(FullGraphUmiPairsIterator other) const;
        bool operator!=(FullGraphUmiPairsIterator other) const;
        std::pair<size_t, size_t> operator*() const;

    private:
        // current_second_ < current_first_
        size_t current_first_;
        size_t current_second_;
        size_t last_;

        void advance();
    };

    class FullGraphUmiPairsIterable {
    public:
        FullGraphUmiPairsIterable(size_t count): count_(count) {}

        FullGraphUmiPairsIterator begin() const;
        FullGraphUmiPairsIterator end() const;

    private:
        size_t count_;
    };


    template <typename ElementType>
    struct Cluster {
    public:
        Cluster(const ElementType& first_, const seqan::Dna5String& center_, const size_t id_, const size_t weight_ = 1)
                : members{first_}, center(center_), id(id_), weight(weight_) {}
        Cluster(const std::unordered_set<ElementType>& members_, const seqan::Dna5String& center_, const size_t weight_, const size_t id_)
                : members(members_), center(center_), id(id_), weight(weight_) {}

        // Assumes merged clusters are never compared.
        bool operator==(const Cluster<ElementType>& other) const { return id == other.id; }

        seqan::Dna5String& GetSequence() const { return center; }

        const std::unordered_set<ElementType> members;
        const seqan::Dna5String center;
        const size_t id;
        const size_t weight;
    };

    template <typename ElementType>
    using ClusterPtr = std::shared_ptr<const Cluster<ElementType>>;

    template <typename ElementType>
    struct ClusterPtrEquals {
        bool operator()(const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second) const {
            return *first == *second;
        }
    };

    template <typename ElementType>
    struct ClusterPtrHash {
        size_t operator()(const ClusterPtr<ElementType>& clusterPtr) const {
            size_t h = clusterPtr->id;
            return h;
        }
    };


    template <typename ElementType>
    using ManyToManyCorrespondenceUmiToCluster = ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>, UmiPtrHash, UmiPtrEquals, ClusterPtrHash<ElementType>, ClusterPtrEquals<ElementType>>;

    template <typename ElementType, typename UmiPairsIterable>
    class Clusterer {
    public:
        static ManyToManyCorrespondenceUmiToCluster<ElementType> cluster(
                const ClusteringMode& mode,
                const std::vector<UmiPtr>& umis,
                const ManyToManyCorrespondenceUmiToCluster<ElementType>& umis_to_clusters,
                const UmiPairsIterable& umi_pairs_iterable);

    private:
        // It's supposed that the merged cluster is not considered in the partition together with either of parents.
        // In particular it inherits id from one of the parents and can be used by ManyToManyCorrespondence, which relies it's a good equality measure and hash value.
        static ClusterPtr<ElementType> mergeClusters(const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second, const size_t id);
        static seqan::Dna5String findNewCenter(std::unordered_set<ElementType>& members);
    };

    template <typename ElementType, typename UmiPairsIterable>
    ManyToManyCorrespondenceUmiToCluster<ElementType> Clusterer<ElementType, UmiPairsIterable>::cluster(
            const ClusteringMode& mode,
            const std::vector<UmiPtr>& umis,
            const ManyToManyCorrespondenceUmiToCluster<ElementType>& umis_to_clusters,
            const UmiPairsIterable& umi_pairs_iterable) {
        ManyToManyCorrespondenceUmiToCluster<ElementType> result(umis_to_clusters);
        // operated on clusters from original umis_to_clusters
        DisjointSets<ClusterPtr<ElementType>> ds;
        for (const auto& cluster : umis_to_clusters.toSet()) {
            ds.addNewSet(cluster);
        }
        std::map<size_t, size_t> dist_distribution;
        for (const auto& umi_pair : umi_pairs_iterable) {
            VERIFY_MSG(umi_pair.first < umis.size() && umi_pair.second < umis.size(), "Invalid umi pair.");
            const UmiPtr& first_umi = umis[umi_pair.first];
            const UmiPtr& second_umi = umis[umi_pair.second];
            for (const auto& first_cluster_original : umis_to_clusters.forth(first_umi)) {
                for (const auto& second_cluster_original : umis_to_clusters.forth(second_umi)) {
                    const auto& first_cluster = result.getTo(ds.findRoot(first_cluster_original));
                    const auto& second_cluster = result.getTo(ds.findRoot(second_cluster_original));
                    if (first_cluster == second_cluster) continue;
                    size_t dist = mode.dist(first_cluster->center, second_cluster->center);
                    dist_distribution[dist] ++;
                    if (dist <= mode.threshold) {
                        // TODO: avoid returning copied umi set by providing access to its begin() and end()
                        const auto& first_cluster_umis = result.back(first_cluster);
                        const auto& second_cluster_umis = result.back(second_cluster);
                        std::unordered_set<UmiPtr> merged_umis(first_cluster_umis);
                        merged_umis.insert(second_cluster_umis.begin(), second_cluster_umis.end());

                        VERIFY_MSG(result.removeTo(first_cluster), "Trying to remove an absent cluster");
                        VERIFY_MSG(result.removeTo(second_cluster), "Trying to remove an absent cluster");

                        VERIFY_MSG(ds.unite(first_cluster_original, second_cluster_original), "Tried to unite two equal sets");
                        const auto merged_cluster = mergeClusters(first_cluster, second_cluster, ds.findRoot(first_cluster_original)->id);

                        result.add(merged_umis, merged_cluster);
                    }
                }
            }
        }

//        INFO("Dist distribution");
//        for (const auto& entry : dist_distribution) {
//            INFO("dist " << entry.first << " : " << entry.second << " times");
//        }

        return result;
    };

    template <typename ElementType, typename UmiPairsIterable>
    ClusterPtr<ElementType> Clusterer<ElementType, UmiPairsIterable>::mergeClusters(const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second, const size_t id) {
        std::unordered_set<ElementType> members(first->members);
        members.insert(second->members.begin(), second->members.end());
        VERIFY_MSG(members.size() == first->members.size() + second->members.size(), "Clusters already intersect.");
        const ClusterPtr<ElementType>& max_cluster = first->members.size() > second->members.size() ? first : second;
        size_t max_size = max_cluster->members.size();
        bool need_new_center = max_size <= 10 || __builtin_clzll(max_size) != __builtin_clzll(members.size());
        const seqan::Dna5String& center = need_new_center ? findNewCenter(members) : max_cluster->center;
        return std::make_shared<Cluster<ElementType>>(members, center, first->weight + second->weight, id);
    }

    template <typename ElementType, typename UmiPairsIterable>
    seqan::Dna5String Clusterer<ElementType, UmiPairsIterable>::findNewCenter(std::unordered_set<ElementType>& members) {
        VERIFY_MSG(members.size() >= 2, "Too few elements to find new center.");
        seqan::Dna5String center = members.begin()->GetSequence();
        std::vector<std::vector<size_t> > cnt(length(center), std::vector<size_t>(4));
        for (auto& member : members) {
            size_t limit = std::min(length(member.GetSequence()), length(center));
            for (size_t i = 0; i < limit; i ++) {
                cnt[i][member.GetSequence()[i]] ++;
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

namespace std {
    template <typename ElementType>
    struct hash<clusterer::Cluster<ElementType>> {
        size_t operator()(const clusterer::Cluster<ElementType>& cluster) const {
            size_t h = cluster.id;
            return h;
        }
    };
}