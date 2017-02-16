#pragma once

#include <seqan/seq_io.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "utils.hpp"
#include "../graph_utils/sparse_graph.hpp"
#include "umi_utils.hpp"
#include "disjoint_sets.hpp"
#include "utils/io.hpp"
#include "../fast_ig_tools/ig_final_alignment.hpp"
#include "../fast_ig_tools/ig_matcher.hpp"

namespace clusterer {

    template <typename ElementType>
    struct Cluster {
    public:
        Cluster(const ElementType& first_, const seqan::Dna5String& center_, const size_t id_, const size_t weight_ = 1)
                : members{first_}, center(center_), id(id_), weight(weight_) {}
        Cluster(const std::unordered_set<ElementType>& members_, const seqan::Dna5String& center_, const size_t id_, const size_t weight_)
                : members(members_), center(center_), id(id_), weight(weight_) {}

        // Assumes merged clusters are never compared.
        bool operator==(const Cluster<ElementType>& other) const { return id == other.id; }

        // TODO: needs recursion once multilevel clustering is needed
        std::unordered_set<Read> GetAllReads() const { return members; }
        size_t size() const { return members.size(); }

        // methods for using Cluster as an ElementType for another Cluster (should be same as Read)
        // TODO: extract this to a trait
        const seqan::Dna5String& GetSequence() const { return center; }

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


    typedef std::function<size_t(const seqan::Dna5String&, const seqan::Dna5String&)> ReadDist;
    typedef std::function<bool(const ClusterPtr<Read>&, const ClusterPtr<Read>&)> ClusterDistChecker;

    struct ClusteringMode {
        ClusteringMode() = delete;

        // return function returning hamming distance if it's <= limit or limit + 1 otherwise
        static ReadDist bounded_hamming_dist(size_t limit);
        // if binary is true, return function returning edit distance if it's <= limit and limit + 1 otherwise
        // else return function returning some number <= limit if edit distance <= limit and + 1 limit otherwise
        static ReadDist bounded_edit_dist(size_t limit, size_t max_indels, bool binary = true);
        static ClusterDistChecker clusters_close_by_center(const ReadDist& read_dist, size_t limit);
        static ClusterDistChecker clusters_close_by_min(const ReadDist& read_dist, size_t limit);
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
                : graph_(graph), vertex_(vertex), current_edge_(current_edge), advances_(0) {}

        GraphUmiPairsIterator operator ++();
        GraphUmiPairsIterator operator ++(int);
        bool operator==(GraphUmiPairsIterator other) const;
        bool operator!=(GraphUmiPairsIterator other) const;
        std::pair<size_t, size_t> operator*() const;

    private:
        const SparseGraphPtr graph_;
        size_t vertex_;
        SparseGraph::EdgesIterator current_edge_;
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
    using ManyToManyCorrespondenceUmiToCluster = ManyToManyCorrespondence<UmiPtr, ClusterPtr<ElementType>, UmiPtrHash, UmiPtrEquals, ClusterPtrHash<ElementType>, ClusterPtrEquals<ElementType>>;

    // ElementType is always Read now. Gradually refactoring to use this information.
    template <typename ElementType>
    class Clusterer {
    public:
        Clusterer(const std::unordered_map<Umi, std::vector<size_t>>& umi_to_reads, const std::unordered_map<Umi, UmiPtr>& umi_ptr_by_umi, const std::vector<Read>& reads);

        template <typename UmiPairsIterable>
        void cluster(const ClusterDistChecker& dist_checker,
                     const std::vector<UmiPtr>& umis,
                     const UmiPairsIterable& umi_pairs_iterable);

        void uniteByCenters();

        void print_umi_split_stats();

        void report_length_differences(const std::string& output_dir);

        void report_non_major_umi_groups_sw(string file_name, string left_graph_file_name, string right_graph_file_name, string chimeras_info_file_name,
                                            string umi_chimeras_info_file_name, size_t num_threads);

        const ManyToManyCorrespondenceUmiToCluster<Read>& getCurrentUmiToCluster() const {
            return current_umi_to_cluster_;
        }

        void write_clusters_and_correspondence(const std::string& base_output_dir, const std::string& postfix, bool save_clusters, bool output,
                                               const std::string& fasta_file_name = "intermediate_repertoire.fasta",
                                               const std::string& rcm_file_name = "intermediate_repertoire.rcm");

    private:
        // It's supposed that the merged cluster is not considered in the partition together with either of parents.
        // In particular it inherits id from one of the parents and can be used by ManyToManyCorrespondence, which relies it's a good equality measure and hash value.
        static ClusterPtr<ElementType> merge_clusters(const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second, const size_t id);
        static seqan::Dna5String findNewCenter(std::unordered_set<ElementType>& members);
        static void print_umi_to_cluster_stats(const ManyToManyCorrespondenceUmiToCluster<ElementType>& umis_to_clusters);

        ManyToManyCorrespondenceUmiToCluster<Read> current_umi_to_cluster_;
        const std::vector<Read> reads_;
    };

    template <typename ElementType>
    Clusterer<ElementType>::Clusterer(const std::unordered_map<Umi, std::vector<size_t>>& umi_to_reads, const std::unordered_map<Umi, UmiPtr>& umi_ptr_by_umi, const std::vector<Read>& reads) :
            reads_(reads) {
        for (auto& entry : umi_to_reads) {
            const auto& umi = umi_ptr_by_umi.at(entry.first);
            for (auto& read_idx : entry.second) {
                const auto& cluster = make_shared<clusterer::Cluster<Read>>(reads_[read_idx], reads_[read_idx].GetSequence(), read_idx);
                current_umi_to_cluster_.add(umi, cluster);
            }
        }
    }

    template <typename ElementType>
    template <typename UmiPairsIterable>
    void Clusterer<ElementType>::cluster(
            const ClusterDistChecker& dist_checker,
            const std::vector<UmiPtr>& umis,
            const UmiPairsIterable& umi_pairs_iterable) {
        ManyToManyCorrespondenceUmiToCluster<Read> result(current_umi_to_cluster_);
        // operated on clusters from original umis_to_clusters
        DisjointSets<ClusterPtr<Read>> ds;
        for (const auto& cluster : current_umi_to_cluster_.toSet()) {
            ds.addNewSet(cluster);
        }
//        std::map<size_t, size_t> dist_distribution;
        for (const auto& umi_pair : umi_pairs_iterable) {
            VERIFY_MSG(umi_pair.first < umis.size() && umi_pair.second < umis.size(), "Invalid umi pair.");
            const UmiPtr& first_umi = umis[umi_pair.first];
            const UmiPtr& second_umi = umis[umi_pair.second];
            for (const auto& first_cluster_original : current_umi_to_cluster_.forth(first_umi)) {
                for (const auto& second_cluster_original : current_umi_to_cluster_.forth(second_umi)) {
                    const auto& first_cluster = result.getTo(ds.findRoot(first_cluster_original));
                    const auto& second_cluster = result.getTo(ds.findRoot(second_cluster_original));
                    if (first_cluster == second_cluster) continue;
                    bool same_cluster = dist_checker(first_cluster, second_cluster);
//                    dist_distribution[dist] ++;
                    if (same_cluster) {
//                    if (dist <= mode.threshold || (dist <= 1.5 * static_cast<double>(mode.threshold) && first_cluster->members.size() == 1 && second_cluster->members.size() == 1)) {
                        // TODO: avoid returning copied umi set by providing access to its begin() and end()
                        const auto& first_cluster_umis = result.back(first_cluster);
                        const auto& second_cluster_umis = result.back(second_cluster);
                        std::unordered_set<UmiPtr, UmiPtrHash, UmiPtrEquals> merged_umis(first_cluster_umis);
                        merged_umis.insert(second_cluster_umis.begin(), second_cluster_umis.end());

                        VERIFY_MSG(result.removeTo(first_cluster), "Trying to remove an absent cluster");
                        VERIFY_MSG(result.removeTo(second_cluster), "Trying to remove an absent cluster");

                        VERIFY_MSG(ds.unite(first_cluster_original, second_cluster_original), "Tried to unite two equal sets");
                        size_t new_id = ds.findRoot(first_cluster_original)->id;
                        const auto merged_cluster = merge_clusters(first_cluster, second_cluster, new_id);

                        result.add(merged_umis, merged_cluster);
                    }
                }
            }
        }

        print_umi_to_cluster_stats(result);

        //        INFO("Dist distribution");
//        for (const auto& entry : dist_distribution) {
//            INFO("dist " << entry.first << " : " << entry.second << " times");
//        }

        current_umi_to_cluster_ = result;
    }

    template <typename ElementType>
    void Clusterer<ElementType>::uniteByCenters() {
        std::unordered_map<seqan::Dna5String, std::vector<ClusterPtr<ElementType>>> center_to_clusters;
        for (const auto& cluster : current_umi_to_cluster_.toSet()) {
            center_to_clusters[cluster->center].push_back(cluster);
        }

        ManyToManyCorrespondenceUmiToCluster<ElementType> result;

        for (const auto& entry : center_to_clusters) {
            const auto& clusters = entry.second;
            std::unordered_set<UmiPtr, UmiPtrHash, UmiPtrEquals> merged_umis(current_umi_to_cluster_.back(clusters[0]));
            ClusterPtr<ElementType> merged_cluster = clusters[0];
            for (size_t i = 1; i < clusters.size(); i ++) {
                const auto& umis = current_umi_to_cluster_.back(clusters[i]);
                merged_umis.insert(umis.begin(), umis.end());
                merged_cluster = merge_clusters(merged_cluster, clusters[i], merged_cluster->id);
            }

            result.add(merged_umis, merged_cluster);
        }

        print_umi_to_cluster_stats(result);

        current_umi_to_cluster_ = result;
    }

    template <typename ElementType>
    void Clusterer<ElementType>::print_umi_split_stats() {
        size_t total_reads = 0;
        size_t secondary_reads = 0;
        for (const auto& umi : current_umi_to_cluster_.fromSet()) {
            const auto& clusters = current_umi_to_cluster_.forth(umi);
            size_t max_cluster = 0;
            size_t umi_size = 0;
            for (const auto& cluster : clusters) {
                umi_size += cluster->weight;
                max_cluster = std::max(max_cluster, cluster->weight);
            }
            total_reads += umi_size;
            secondary_reads += umi_size - max_cluster;
        }
        INFO("Total " << total_reads << " reads in the dataset. " << secondary_reads << " of them don't belong to the main group of their barcode.");
    }

    template <typename ElementType>
    void Clusterer<ElementType>::report_length_differences(const std::string& output_dir) {
        namespace fs = boost::filesystem;
        if (output_dir != ".") {
            fs::create_directory(output_dir);
            INFO("Created new");
        }

        std::ofstream file_with_seqs(fs::path(output_dir).append("clusters_with_huge_len_span.txt").string());

        map<size_t, size_t> len_span_to_count;
        map<size_t, map<size_t, size_t>> len_span_to_sizes;
        for (const auto& cluster : current_umi_to_cluster_.toSet()) {
            size_t min_len = numeric_limits<size_t>::max();
            size_t max_len = 0;
            for (const auto& read : cluster->members) {
                size_t len = length(read.GetSequence());
                min_len = min(min_len, len);
                max_len = max(max_len, len);
            }
            len_span_to_count[max_len - min_len] ++;
            len_span_to_sizes[max_len - min_len][cluster->members.size()] ++;

            if (max_len - min_len > 10) {
                file_with_seqs << (max_len - min_len) << "\n";
                for (const auto& read : cluster->members) {
                    file_with_seqs << ">" << read.GetReadId() << "\n" << read.GetSequence() << "\n";
                }
            }
        }
        INFO("Number of clusters by difference between shortest and longest member.")
        for (const auto& entry : len_span_to_count) {
            stringstream ss;
            ss << entry.first << " -> " << entry.second << ". ";
            for (const auto& size_entry : len_span_to_sizes[entry.first]) {
                ss << " " << size_entry.first << ": " << size_entry.second;
            }
            INFO(ss.str());
        }
    }

    template <typename ElementType>
    void Clusterer<ElementType>::report_non_major_umi_groups_sw(std::string file_name,
                                                   std::string left_graph_file_name,
                                                   std::string right_graph_file_name,
                                                   std::string chimeras_info_file_name,
                                                   std::string umi_chimeras_info_file_name,
                                                   size_t num_threads) {
        // TODO: 1) 10 is ok to be close, but too few to be different; 20 should probably do.
        // TODO: 2) Too slow. First compute candidates inside UMI, even quadratically, then use repr k-mers for another half.
        // TODO:    Another option is to check if actually 50-mer strategy is precise enough.
        std::vector<seqan::Dna5String> all_left_halves;
        std::vector<seqan::Dna5String> all_right_halves;
        const size_t IG_LEN = 350;
        std::unordered_map<seqan::Dna5String, size_t> cluster_to_idx;
        for (const auto& cluster : current_umi_to_cluster_.toSet()) {
            const auto& sequence = cluster->GetSequence();
            cluster_to_idx[cluster->GetSequence()] = all_left_halves.size();
            all_left_halves.push_back(seqan::prefix(sequence, IG_LEN / 2));
            all_right_halves.push_back(seqan::suffix(sequence, IG_LEN / 2));
        }

        const size_t tau = 10;
        const size_t max_indels = 3;
        const size_t strategy = 2;
        const size_t k = IG_LEN / (tau + strategy) / 2;
        Graph left_graph;
        Graph right_graph;
        omp_set_num_threads(static_cast<int>(num_threads));
        for (size_t i = 0; i < 2; i ++) {
            INFO("constructing graph " << i);
            const std::vector<seqan::Dna5String>& all_halves = (i == 0) ? all_left_halves : all_right_halves;
            auto kmer2reads = kmerIndexConstruction(all_halves, k);
            size_t num_of_dist_computations = 0;
            Graph& graph = (i == 0) ? left_graph : right_graph;
            graph = tauDistGraph(all_halves, kmer2reads, ClusteringMode::bounded_edit_dist(tau, max_indels),
                                 static_cast<unsigned>(tau), static_cast<unsigned>(k), static_cast<unsigned>(strategy),
                                 num_of_dist_computations);
            INFO("graph constructed: " << i << ", dist computed " << num_of_dist_computations << " times.");
            write_metis_graph(graph, i == 0 ? left_graph_file_name : right_graph_file_name);
            INFO("graph written: " << i);
        }

        for (size_t i = 0; i < left_graph.size(); i ++) {
            std::sort(left_graph[i].begin(), left_graph[i].end());
            std::sort(right_graph[i].begin(), right_graph[i].end());
        }

        size_t total_minor_clusters = 0;
        size_t found_somewhere = 0;
        size_t found_within_umi = 0;
        size_t chimeric_clusters = 0;
        size_t chimeric_reads = 0;
        size_t found_half_only = 0;
        std::map<size_t, size_t> chimera_size_to_count;

        std::ofstream out_file(file_name);
        std::ofstream chimeras_file(chimeras_info_file_name);
        std::ofstream umi_chimeras_file(umi_chimeras_info_file_name);
        auto result = current_umi_to_cluster_;
        const auto& umis = current_umi_to_cluster_.fromSet();
        for (const auto& umi: umis) {
            const auto& clusters = current_umi_to_cluster_.forth(umi);
            size_t max_size = 0;
            seqan::Dna5String max_consensus;
            for (const auto& cluster : clusters) {
                if (cluster->weight > max_size) {
                    max_size = cluster->weight;
                    max_consensus = cluster->GetSequence();
                }
            }
            // counting major UMI clusters only
            std::unordered_set<size_t> major_umi_clusters_idcs;
            for (const auto& cluster : clusters) {
                if (cluster->weight == max_size) {
                    major_umi_clusters_idcs.insert(cluster_to_idx[cluster->GetSequence()]);
                }
            }

            bool was_max = false;
            for (const auto& cluster : clusters) {
                if (cluster->weight == max_size) {
                    if (was_max) {
                        out_file << cluster->weight << "\t" << max_size << "\tp" << std::endl;
                        if (max_size > 1) {
                            INFO("Max component parity: " << cluster->weight);
                        }
                    }
                    was_max = true;
                } else {
                    out_file << cluster->weight << "\t" << max_size << std::endl;
                }

                if (max_size < 3) continue;
                if (cluster->weight == max_size) continue;

                size_t left_in_umi = 0;
                const auto& left_adjacent = left_graph[cluster_to_idx[cluster->GetSequence()]];
                size_t left_in_all = left_adjacent.size();
                for (size_t i = 0; i < left_adjacent.size(); i ++) {
                    if (major_umi_clusters_idcs.count(left_adjacent[i].first)) {
                        left_in_umi ++;
                    }
                }
                size_t right_in_umi = 0;
                const auto& right_adjacent = right_graph[cluster_to_idx[cluster->GetSequence()]];
                size_t right_in_all = right_adjacent.size();
                for (size_t i = 0; i < right_adjacent.size(); i ++) {
                    if (major_umi_clusters_idcs.count(right_adjacent[i].first)) {
                        right_in_umi ++;
                    }
                }

                total_minor_clusters ++;
                if (left_in_all > 0 && right_in_all > 0) {
                    found_somewhere ++;
                }
                const auto sw_dist = ClusteringMode::bounded_edit_dist(40, 15, false);
                if (left_in_umi > 0 && right_in_umi > 0) {
                    found_within_umi ++;
                    umi_chimeras_file << "size: " << cluster->weight << " (max = " << max_size << ")" << std::endl;
                    umi_chimeras_file << "chimera?: " << cluster->GetSequence() << std::endl;
                    umi_chimeras_file << "max consensus: " << max_consensus << std::endl;
                    umi_chimeras_file << left_in_all << " " << right_in_all << " " << left_in_umi << " " << right_in_umi << std::endl;
                    umi_chimeras_file << "umi clusters with dists:" << std::endl;
                    for (const auto& c : clusters) {
                        umi_chimeras_file << "cluster size: " << c->weight << ", ";
                        umi_chimeras_file << "left dist: " << sw_dist(seqan::prefix(cluster->GetSequence(), IG_LEN / 2), seqan::prefix(c->GetSequence(), IG_LEN / 2)) << ", ";
                        umi_chimeras_file << "right dist: " << sw_dist(seqan::suffix(cluster->GetSequence(), IG_LEN / 2), seqan::suffix(c->GetSequence(), IG_LEN / 2)) << std::endl;
                        umi_chimeras_file << c->GetSequence() << std::endl;
                    }
                } else if ((left_in_umi > 0 || right_in_umi > 0) && left_in_all > 0 && right_in_all > 0) {
                    bool is_chimera = true;
//                    for (const auto& c : clusters) {
//                        // This is not a chimera, but a highly amplified sequence
//                        if (sw_dist(cluster->GetSequence(), c->GetSequence()) <= 25) {
//                            is_chimera = false;
//                            break;
//                        }
//                    }
                    if (!is_chimera) continue;

                    chimeric_reads += cluster->members.size();
                    chimera_size_to_count[cluster->members.size()] ++;
                    result.removeTo(cluster);
                    chimeric_clusters ++;
                    chimeras_file << "size: " << cluster->weight << " (max = " << max_size << ")" << std::endl;
                    chimeras_file << "chimera: " << cluster->GetSequence() << std::endl;
                    chimeras_file << "max consensus: " << max_consensus << std::endl;
                    chimeras_file << left_in_all << " " << right_in_all << " " << left_in_umi << " " << right_in_umi << std::endl;
                    chimeras_file << "umi clusters with dists:" << std::endl;
                    for (const auto& c : clusters) {
                        chimeras_file << "cluster size: " << c->weight << ", ";
                        chimeras_file << "left dist: " << sw_dist(seqan::prefix(cluster->GetSequence(), IG_LEN / 2), seqan::prefix(c->GetSequence(), IG_LEN / 2)) << ", ";
                        chimeras_file << "right dist: " << sw_dist(seqan::suffix(cluster->GetSequence(), IG_LEN / 2), seqan::suffix(c->GetSequence(), IG_LEN / 2)) << std::endl;
                        chimeras_file << c->GetSequence() << std::endl;
                    }
                } else if ((left_in_all > 0) != (right_in_all > 0)) {
                    found_half_only ++;
                }
            }
        }

        INFO("Total clusters, which are minorities in the umi: " << total_minor_clusters);
        INFO("From them could be found elsewhere by halves (maybe in the same source): " << found_somewhere);
        INFO("Ones that could be found within the same UMI by halves: " << found_within_umi);
        INFO("Clusters with chimeras - ones with a half within UMI and a half somewhere else outside: " << chimeric_clusters);
        INFO("Ones with only one half found at all: " << found_half_only);
        INFO("Total chimeric reads: " << chimeric_reads);
        INFO("Chimera size distribution (size -> count):");
        for (const auto& entry : chimera_size_to_count) {
            INFO(entry.first << " -> " << entry.second);
        }

        current_umi_to_cluster_ = result;
    }

    template <typename ElementType>
    void Clusterer<ElementType>::write_clusters_and_correspondence(const std::string& base_output_dir, const std::string& postfix, bool save_clusters, bool output,
                                                                   const std::string& fasta_file_name, const std::string& rcm_file_name) {
        if (!output) return;

        namespace fs = boost::filesystem;

        auto abs_base_output_dir = fs::absolute(fs::path(base_output_dir)).normalize();
        while (abs_base_output_dir.filename() == ".") {
            abs_base_output_dir = abs_base_output_dir.parent_path();
        }
        const auto output_dir = abs_base_output_dir.parent_path() / (abs_base_output_dir.filename().string() + postfix);

        std::vector<seqan::CharString> repertoire_ids;
        std::vector<seqan::Dna5String> repertoire_reads;
        std::unordered_map<size_t, size_t> read_id_to_cluster_id;
        {
            size_t cluster_id = 0;
            for (const auto& cluster : current_umi_to_cluster_.toSet()) {
                repertoire_ids.emplace_back("intermediate_cluster___" + to_string(cluster_id) + "___size___" + to_string(cluster->size()));
                repertoire_reads.push_back(cluster->center);
                VERIFY(cluster->GetAllReads().size() > 0);
                for (const auto& read : cluster->GetAllReads()) {
                    read_id_to_cluster_id[read.GetId()] = cluster_id;
                }
                cluster_id++;
            }
        }
        VERIFY(read_id_to_cluster_id.size() <= reads_.size());

        INFO("Using output directory " << output_dir);
        if (output_dir != ".") {
            fs::create_directory(output_dir);
            INFO("Created new");
        }
        if (!fasta_file_name.empty()) {
            write_seqan_records(fs::path(output_dir).append(fasta_file_name), repertoire_ids, repertoire_reads);
        }

        std::ofstream read_to_cluster_ofs(fs::path(output_dir).append(rcm_file_name).string());
        for (const auto& read : reads_) {
            read_to_cluster_ofs << read.GetReadId();
            if (read_id_to_cluster_id.count(read.GetId())) {
                read_to_cluster_ofs << "\t" << read_id_to_cluster_id[read.GetId()];
            }
            read_to_cluster_ofs << "\n";
        }
        read_to_cluster_ofs.close();

        if (!save_clusters) return;

        const auto clusters_path = fs::path(output_dir).append("clusters_by_umis");
        create_new_directory(clusters_path);
        for (const auto& umi : current_umi_to_cluster_.fromSet()) {
            auto umi_path = clusters_path;
            umi_path.append(seqan_string_to_string(umi->GetString()));
            fs::create_directory(umi_path);
            const auto& cluster_set = current_umi_to_cluster_.forth(umi);
            std::vector<ClusterPtr<ElementType>> clusters(cluster_set.begin(), cluster_set.end());
            std::sort(clusters.begin(), clusters.end(), [](const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second) {
                if (first->weight != second->weight) {
                    return first->weight > second->weight;
                }
                return first->id < second->id;
            });
            for (const auto& cluster : clusters) {
                std::vector<seqan::CharString> cluster_ids;
                std::vector<seqan::Dna5String> cluster_reads;
                for (const auto& read : cluster->members) {
                    cluster_ids.push_back(read.GetReadId());
                    cluster_reads.push_back(read.GetSequence());
                }
                const auto cluster_path = fs::path(umi_path).append("cluster_" + std::to_string(cluster->id) + "_size_" + std::to_string(cluster->weight) + ".fasta");
                write_seqan_records(cluster_path, cluster_ids, cluster_reads);
            }
        }
    }

    template <typename ElementType>
    void Clusterer<ElementType>::print_umi_to_cluster_stats(const ManyToManyCorrespondenceUmiToCluster<ElementType>& umis_to_clusters) {
        map<size_t, size_t> clusters_per_umi;
        for (const auto& umi : umis_to_clusters.fromSet()) {
            clusters_per_umi[umis_to_clusters.forth(umi).size()] ++;
        }
        INFO("Distribution of number of clusters per UMI: size -> count");
        for (const auto& entry : clusters_per_umi) {
            INFO(entry.first << "\t" << entry.second);
        }

        std::map<size_t, size_t> umis_per_cluster;
        for (const auto& cluster : umis_to_clusters.toSet()) {
            umis_per_cluster[umis_to_clusters.back(cluster).size()] ++;
        }
        INFO("Distribution of number of UMIs per cluster: size -> count");
        for ( const auto& entry : umis_per_cluster) {
            INFO(entry.first << "\t" << entry.second);
        }
    };

    template <typename ElementType>
    ClusterPtr<ElementType> Clusterer<ElementType>::merge_clusters(const ClusterPtr<ElementType>& first, const ClusterPtr<ElementType>& second, const size_t id) {
        std::unordered_set<ElementType> members(first->members);
        members.insert(second->members.begin(), second->members.end());
        VERIFY_MSG(members.size() == first->members.size() + second->members.size(), "Clusters already intersect.");
        const ClusterPtr<ElementType>& max_cluster = first->members.size() > second->members.size() ? first : second;
        size_t max_size = max_cluster->members.size();
        bool need_new_center = max_size <= 20 || __builtin_clzll(max_size) != __builtin_clzll(members.size());
        const seqan::Dna5String& center = need_new_center ? findNewCenter(members) : max_cluster->center;
        return std::make_shared<Cluster<ElementType>>(members, center, id, first->weight + second->weight);
    }

    template <typename ElementType>
    seqan::Dna5String Clusterer<ElementType>::findNewCenter(std::unordered_set<ElementType>& members) {
        VERIFY_MSG(members.size() >= 2, "Too few elements to find new center.");
        std::vector<seqan::Dna5String> reads(members.size());
        std::vector<size_t> indices(members.size());
        size_t current = 0;
        for (Read member : members) {
            reads[current] = member.GetSequence();
            indices[current] = current;
            current ++;
        }
        return consensus_hamming_limited_coverage(reads, indices, std::vector<size_t>(reads.size(), 1), 1);
        /*seqan::Dna5String center = seqan::Dna5String(members.begin()->GetSequence());
        std::vector<std::vector<size_t>> cnt(length(center), std::vector<size_t>(4));
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
        return center;*/
    }

    template <typename ElementType>
    size_t count_reads_with_corrected_umi(const ManyToManyCorrespondenceUmiToCluster<ElementType>& umi_to_clusters) {
        const auto& clusters = umi_to_clusters.toSet();
        size_t result = 0;
//        size_t printed = 0;
        for (const auto& cluster : clusters) {
            const auto& members = cluster->members;
            std::vector<seqan::CharString> read_ids(members.size());
            std::transform(members.begin(), members.end(), read_ids.begin(), [](Read read) -> seqan::CharString { return read.GetReadId(); });
            std::vector<seqan::Dna5String> umis;
            extract_barcodes_from_read_ids(read_ids, umis);
            std::unordered_map<seqan::Dna5String, size_t> umi_to_cnt;
            for (const seqan::Dna5String& umi : umis) {
                umi_to_cnt[umi] ++;
            }
            size_t max_cnt = 0;
            seqan::Dna5String max_umi;
            for (const auto& entry : umi_to_cnt) {
                if (entry.second > max_cnt) {
                    max_cnt = entry.second;
                    max_umi = entry.first;
                }
            }
            bool was_max = false;
            for (const auto& entry : umi_to_cnt) {
                if (entry.second == max_cnt) {
                    if (was_max) {
                        result += entry.second;
                        TRACE("Tie for count " << max_cnt);
                    }
                    was_max = true;
                } else {
                    result += entry.second;
                }
            }

//            for (const auto& read: members) {
//                if (seqan_string_to_string(read.GetReadId()).find(seqan_string_to_string(max_umi)) == std::string::npos) {
//                    std::cout << ">" << read.GetReadId() << std::endl << read.GetSequence() << std::endl;
//                    printed ++;
//                }
//            }
        }
//        if (result != printed) {
//            ERROR("Returning " << result << ", but printed " << printed);
//        }
        return result;
    }

    template <typename ElementType>
    void report_non_major_umi_groups(const ManyToManyCorrespondenceUmiToCluster<ElementType> umi_to_clusters, std::string file_name) {
        std::unordered_map<std::string, size_t> all_left_kmers;
        std::unordered_map<std::string, size_t> all_right_kmers;
        const size_t IG_LEN = 350;
        const size_t K = IG_LEN / 6;
        const size_t MARGIN = (IG_LEN / 2 - K) / 2;
        for (const auto& cluster : umi_to_clusters.toSet()) {
//            for (const Read& read : cluster->members) {
//                std::string sequence = seqan_string_to_string(read.GetSequence());
//                for (size_t shift = 0; shift <= MARGIN * 2; shift ++) {
//                    all_left_kmers[sequence.substr(shift, K)] ++;
//                    all_right_kmers[sequence.substr(sequence.length() - K - shift, K)] ++;
//                }
//            }
            std::string sequence = seqan_string_to_string(cluster->GetSequence());
            for (size_t shift = 0; shift <= MARGIN * 2; shift ++) {
                all_left_kmers[sequence.substr(shift, K)] ++;
                all_right_kmers[sequence.substr(sequence.length() - K - shift, K)] ++;
            }
        }

        size_t total_singles_minor = 0;
        size_t found_somewhere = 0;
        size_t found_within_umi = 0;
        size_t found_partially = 0;
//        size_t found_partially_left = 0;
//        size_t found_partially_right = 0;

        std::ofstream out_file(file_name);
        for (const auto& umi: umi_to_clusters.fromSet()) {
            std::unordered_map<std::string, size_t> umi_left_kmers;
            std::unordered_map<std::string, size_t> umi_right_kmers;
            for (const auto& cluster : umi_to_clusters.forth(umi)) {
//                for (const Read& read : cluster->members) {
//                    std::string sequence = seqan_string_to_string(read.GetSequence());
//                    for (size_t shift = 0; shift <= MARGIN * 2; shift ++) {
//                        umi_left_kmers[sequence.substr(shift, K)] ++;
//                        umi_right_kmers[sequence.substr(sequence.length() - K - shift, K)] ++;
//                    }
//                }
                std::string sequence = seqan_string_to_string(cluster->GetSequence());
                for (size_t shift = 0; shift <= MARGIN * 2; shift ++) {
                    all_left_kmers[sequence.substr(shift, K)] ++;
                    all_right_kmers[sequence.substr(sequence.length() - K - shift, K)] ++;
                }
            }

            size_t max_size = 0;
            seqan::Dna5String max_consensus;
            const auto& clusters = umi_to_clusters.forth(umi);
            for (const auto& cluster : clusters) {
                if (cluster->weight > max_size) {
                    max_size = cluster->weight;
                    max_consensus = cluster->GetSequence();
                }
            }
            bool was_max = false;
            for (const auto& cluster : clusters) {
                if (cluster->weight == max_size) {
                    if (was_max) {
                        out_file << cluster->weight << "\t" << max_size << "\tp" << std::endl;
                        INFO("Max component parity: " << cluster->weight);
                    }
                    was_max = true;
                } else {
                    out_file << cluster->weight << "\t" << max_size << std::endl;
                }
                if (max_size > 1 && cluster->weight < max_size /*> 2*/) {
                    total_singles_minor ++;
                    std::string sequence = seqan_string_to_string(cluster->members.begin()->GetSequence());
                    std::string left_kmer = sequence.substr(MARGIN, K);
                    std::string right_kmer = sequence.substr(sequence.length() - K - MARGIN, K);
                    if (all_left_kmers[left_kmer] > 1 && all_right_kmers[right_kmer] > 1) {
                        found_somewhere ++;
                    }
                    if (umi_left_kmers[left_kmer] > 1 && umi_right_kmers[right_kmer] > 1) {
                        found_within_umi ++;
                    } else if ((umi_left_kmers[left_kmer] > 1 || umi_right_kmers[right_kmer] > 1) &&
                            (all_left_kmers[left_kmer] > 1 || all_right_kmers[right_kmer] > 1)) {
                        found_partially ++;
                        INFO("size: " << cluster->weight);
                        INFO("chimera?: " << sequence);
                        INFO("consensus: " << max_consensus);
                        INFO(all_left_kmers[left_kmer] << " " << all_right_kmers[right_kmer] << " " << umi_left_kmers[left_kmer] << " " << umi_right_kmers[right_kmer])
                    }

                }
            }
        }

        INFO("Total singletons, which are minorities in the umi: " << total_singles_minor);
        INFO("From them could be found elsewhere by halves (maybe in the same source): " << found_somewhere);
        INFO("Ones that could be found within the same UMI by halves: " << found_within_umi);
        INFO("Ones with a half within UMI and a half somewhere else outside: " << found_partially);
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

template <typename T>
T abs_diff(T first, T second) {
    return first > second ? first - second : second - first;
}
