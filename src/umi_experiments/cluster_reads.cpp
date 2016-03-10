#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <seqan/seq_io.h>
#include "utils.hpp"
#include "umi_utils.hpp"
#include "../graph_utils/graph_io.hpp"
#include "disjoint_sets.hpp"
#include "clusterer.hpp"

struct Params {
    std::string reads_path;
    std::string umi_uncompressed_path;
    std::string umi_compressed_path;
    std::string umi_graph_path;
};

bool read_args(int argc, char **argv, Params& params) {
    namespace po = boost::program_options;
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("reads,r", po::value<std::string>(&params.reads_path)->required(), "input file with reads")
            ("umi_uncompressed,u", po::value<std::string>(&params.umi_uncompressed_path)->required(), "file with UMI records extracted (not compressed)")
            ("umi_compressed,c", po::value<std::string>(&params.umi_compressed_path)->required(), "file with UMI records extracted (compressed)")
            ("graph,g", po::value<std::string>(&params.umi_graph_path)->required(), "file with UMI graph")
//            ("output,o", po::value<std::string>(&output_dir)->default_value(""), "output directory")
//            ("threads,t", po::value<unsigned>(&thread_count)->default_value(1), "number of running threads")
            ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
    if (vm.count("help") || argc == 1) {
        std::cout << cmdl_options << std::endl;
        return false;
    }
    po::notify(vm);
    return true;
}

class Clusterer {
public:
    class Cluster {
    public:
        Cluster() {}

        Cluster(size_t id, const seqan::Dna5String& first_member, size_t first_member_index)
                : id_(id), members_{first_member_index}, center_(first_member) {}

        Cluster(size_t id, const std::vector<size_t>& members, const seqan::Dna5String& dummy_center, const std::vector<seqan::Dna5String> &reads)
                : id_(id), members_(members), center_(dummy_center) {
            reevalueteCenter(reads);
        }

        bool shouldContain(const seqan::Dna5String &read, const clusterer::ClusteringMode &clusteringMode) const {
            return clusteringMode.dist(read, center_) <= clusteringMode.threshold;
        }

        void add(const std::vector<seqan::Dna5String> &reads, size_t index) {
            members_.push_back(index);
            unsigned long size = members_.size();
            if (size < 10 || ((size & (size - 1)) == 0)) {
                reevalueteCenter(reads);
            }
        }

        static Cluster merge(const Cluster& first, const Cluster& second, const std::vector<seqan::Dna5String> &reads) {
            std::vector<size_t> members;
            members.insert(members.end(), first.members_.begin(), first.members_.end());
            members.insert(members.end(), second.members_.begin(), second.members_.end());
            return Cluster(first.id_, members, first.center_, reads);
        }

        // Larger clusters go first
        struct SizeComparator {
            bool operator() (const Cluster& left, const Cluster& right) const {
                if (left.members_.size() > right.members_.size()) {
                    return true;
                }
                return left.members_.size() == right.members_.size() && left.id_ > right.id_;
            }
        };

    private:
        size_t id_;
        std::vector<size_t> members_;
        seqan::Dna5String center_;

        void reevalueteCenter(const std::vector<seqan::Dna5String> &reads) {
            std::vector<std::vector<size_t> > cnt(length(center_), std::vector<size_t>(4));
            for (auto index : members_) {
                size_t limit = std::min(length(reads[index]), length(center_));
                for (size_t i = 0; i < limit; i ++) {
                    cnt[i][reads[index][i]] ++;
                }
            }
            for (size_t i = 0; i < length(center_); i ++) {
                for (size_t j = 0; j < cnt[i].size(); j ++) {
                    if (cnt[i][j] > cnt[i][center_[i]]) {
                        center_[i] = j;
                    }
                }
            }
        }
    };

    Clusterer(const clusterer::ClusteringMode& mode) : mode_(mode) {}

    std::vector<Cluster> cluster(const std::vector<size_t>& indices, const std::vector<seqan::Dna5String>& reads) {
        std::set<Cluster, Cluster::SizeComparator> clusters;
        for (auto index : indices) {
            bool found = false;
            for (const auto& cluster : clusters) {
                if (cluster.shouldContain(reads[index], mode_)) {
                    auto new_cluster = cluster;
                    clusters.erase(cluster);
                    new_cluster.add(reads, index);
                    clusters.insert(new_cluster);
                    found = true;
                    break;
                }
            }
            if (!found) {
                clusters.insert(Cluster(index, reads[index], index));
            }
        }
        return std::vector<Cluster>(clusters.begin(), clusters.end());
    }

private:
    const clusterer::ClusteringMode mode_;
};

size_t cluster_inside_umi(const std::vector<seqan::Dna5String>& input_reads,
                          const std::unordered_map <Umi, std::vector<size_t>>& umi_to_reads,
                          std::unordered_map <Umi, std::vector<Clusterer::Cluster>>& hamm_clusters_by_umi) {
    size_t total_clusters = 0;
    for (auto& entry : umi_to_reads) {
        const auto& umi = entry.first;
        const auto& umi_reads = entry.second;
        hamm_clusters_by_umi[umi] = Clusterer(clusterer::ClusteringMode::hamming).cluster(umi_reads, input_reads);
        total_clusters += hamm_clusters_by_umi[umi].size();
    }
    return total_clusters;
}

void unite_clusters_for_adjacent_umis(const std::vector<seqan::Dna5String>& input_reads, const SparseGraphPtr umi_graph,
                                        const std::vector<seqan::Dna5String>& umis,
                                        const std::unordered_map<Umi, std::vector<Clusterer::Cluster>>& hamm_clusters_by_umi,
                                        std::vector<Clusterer::Cluster>& hamm_clusters) {
    // TODO: share Umi
    typedef std::pair<Umi, size_t> ClusterId;
    struct ClusterIdHash {
        size_t operator()(const ClusterId& id) const {
            return std::hash<Umi>()(id.first) ^ std::hash<size_t>()(id.second);
        }
    };
    DisjointSets<ClusterId, ClusterIdHash> cluster_origins;
    std::unordered_map<ClusterId, Clusterer::Cluster, ClusterIdHash> clusters;
    for (auto& entry : hamm_clusters_by_umi) {
        auto& umi = entry.first;
        auto& umi_clusters = entry.second;
        for (size_t i = 0; i < umi_clusters.size(); i ++) {
            auto cluster_id = std::make_pair(umi, i);
            cluster_origins.addNewSet(cluster_id);
            clusters[cluster_id] = umi_clusters[i];
        }
    }

    for (size_t v = 0; v < umi_graph->N(); v ++) {
        const auto v_umi = Umi(umis[v]);
        const auto& v_clusters = hamm_clusters_by_umi.find(v_umi)->second;
        for (size_t u : umi_graph->VertexEdges(v)) {
            const auto u_umi = Umi(umis[u]);
            const auto& u_clusters = hamm_clusters_by_umi.find(u_umi)->second;
            for (size_t v_cluster_idx = 0; v_cluster_idx < v_clusters.size(); v_cluster_idx ++) {
                const auto v_cluster_id = std::make_pair(v_umi, v_cluster_idx);
                for (size_t u_cluster_idx = 0; u_cluster_idx < u_clusters.size(); u_cluster_idx ++) {
                    const auto u_cluster_id = std::make_pair(u_umi, u_cluster_idx);
                    if (clusterer::ClusteringMode::hamming.dist(umis[u], umis[v]) > clusterer::ClusteringMode::hamming.threshold) {
                        continue;
                    }
                    if (!cluster_origins.unite(v_cluster_id, u_cluster_id)) {
                        continue;
                    }
                    const auto& merged = Clusterer::Cluster::merge(v_clusters[v_cluster_idx], u_clusters[u_cluster_idx], input_reads);
                    clusters.erase(v_cluster_id);
                    clusters.erase(u_cluster_id);
                    clusters[cluster_origins.findRoot(v_cluster_id)] = merged;
                }
            }
        }
    }

    for (auto& entry : clusters) {
        hamm_clusters.push_back(entry.second);
    }
}

struct Input {
    Input(vector<seqan::CharString>& input_ids_, std::vector<seqan::Dna5String>& input_reads_,
          vector<seqan::CharString>& umi_ids_, std::vector<seqan::Dna5String>& umis_,
          std::vector<seqan::Dna5String>& compressed_umis_, SparseGraphPtr& umi_graph_)
            : input_ids(input_ids_), input_reads(input_reads_), umi_ids(umi_ids_), umis(umis_),
              compressed_umis(compressed_umis_), umi_graph(umi_graph_) {}

    vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    vector<seqan::CharString> umi_ids;
    std::vector<seqan::Dna5String> umis;
    std::vector<seqan::Dna5String> compressed_umis;
    SparseGraphPtr umi_graph;
};

Input read_everything(const Params& params) {
    vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    INFO("Reading reads");
    seqan::SeqFileIn reads_file(params.reads_path.c_str());
    readRecords(input_ids, input_reads, reads_file);
    INFO(input_ids.size() << " records read");

    vector<seqan::CharString> umi_ids;
    std::vector<seqan::Dna5String> umis;
    INFO("Reading UMI records");
    seqan::SeqFileIn umi_file(params.umi_uncompressed_path.c_str());
    readRecords(umi_ids, umis, umi_file);
    INFO(umi_ids.size() << " UMIs read");

    std::vector<seqan::Dna5String> compressed_umis;
    INFO("Reading compressed UMIs");
    std::vector<seqan::CharString> compressed_umi_ids;
    seqan::SeqFileIn umi_compressed_file(params.umi_compressed_path.c_str());
    readRecords(compressed_umi_ids, compressed_umis, umi_compressed_file);
    INFO(compressed_umis.size() << " compressed UMIs read");

    SparseGraphPtr umi_graph;
    INFO("Reading UMI graph");
    umi_graph = GraphReader(params.umi_graph_path).CreateGraph();
    INFO("Read graph with " << umi_graph->N() << " vertices and " << umi_graph->NZ() << " edges");

    return Input(input_ids, input_reads, umi_ids, umis, compressed_umis, umi_graph);
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    const auto input = read_everything(params);

    // needs uncompressed umis
    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(input.umis, umi_to_reads);


    INFO("Employing new sctructure");
    // TODO: get rid of plain Dna5Strings in favor of shared_ptrs to them
    std::unordered_map<Umi, UmiPtr> umi_ptr_by_umi;
    for (auto& umi_read : input.compressed_umis) {
        umi_ptr_by_umi[Umi(umi_read)] = std::make_shared<Umi>(umi_read);
    }

    std::vector<UmiPtr> compressed_umi_ptrs;
    for (auto& entry : umi_ptr_by_umi) {
        compressed_umi_ptrs.push_back(entry.second);
    }
    clusterer::ManyToManyCorrespondenceUmiToCluster<Read> initial_umis_to_clusters;
    for (auto& entry : umi_to_reads) {
        auto& umi = umi_ptr_by_umi[entry.first];
        for (auto& read_idx : entry.second) {
            const auto& read = Read(input.input_reads[read_idx], read_idx);
            const auto& cluster = std::make_shared<clusterer::Cluster<Read>>(read, input.input_reads[read_idx], read_idx);
            initial_umis_to_clusters.add(umi, cluster);
        }
    }
    INFO("Clustering reads by hamming within single UMIs with threshold " << clusterer::ClusteringMode::hamming.threshold);
    const auto& umi_to_clusters_hamm_inside_umi = clusterer::Clusterer<Read, clusterer::ReflexiveUmiPairsIterable>::cluster(
            clusterer::ClusteringMode::hamming, compressed_umi_ptrs, initial_umis_to_clusters,
            clusterer::ReflexiveUmiPairsIterable(compressed_umi_ptrs.size()));
    INFO(umi_to_clusters_hamm_inside_umi.toSize() << " clusters found");


    INFO("Old");
    INFO("Clustering reads by hamming within single UMIs with threshold " << clusterer::ClusteringMode::hamming.threshold);
    std::unordered_map<Umi, std::vector<Clusterer::Cluster>> hamm_clusters_by_umi;
    {
        size_t total_clusters = cluster_inside_umi(input.input_reads, umi_to_reads, hamm_clusters_by_umi);
        INFO(total_clusters << " clusters found");
    }

    // unite groups of close by hamming reads for adjacent UMIs
    INFO("Uniting read clusters for adjacent UMIs");
    std::vector<Clusterer::Cluster> hamm_clusters;
    unite_clusters_for_adjacent_umis(input.input_reads, input.umi_graph, input.compressed_umis, hamm_clusters_by_umi, hamm_clusters);
    INFO(hamm_clusters.size() << " clusters found");

    

    // probably proceed with edit distance

    // unite close reads with different UMIs: graph is needed anyway; then either metis clustering, or continue custom techniques

    return 0;
}
