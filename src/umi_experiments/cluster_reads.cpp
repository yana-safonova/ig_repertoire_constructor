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

namespace {
    struct Params {
        std::string reads_path;
        std::string umi_uncompressed_path;
        std::string umi_compressed_path;
        std::string umi_graph_path;
        std::string output_dir;
        bool save_clusters;
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
                ("output,o", po::value<std::string>(&params.output_dir)->default_value(""), "output directory path")
                ("save-clusters,s", po::value<bool>(&params.save_clusters)->default_value(false), "save clusters by UMI")
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

    struct Input {
        Input(vector<seqan::CharString>& input_ids_, std::vector<seqan::Dna5String>& input_reads_,
              vector<seqan::CharString>& umi_ids_, std::vector<seqan::Dna5String>& umis_,
              std::vector<seqan::Dna5String>& compressed_umis_, SparseGraphPtr& umi_graph_)
                : input_ids(input_ids_), input_reads(input_reads_), umi_ids(umi_ids_), umis(umis_),
                  compressed_umis(compressed_umis_), umi_graph(umi_graph_) {}

        std::vector<seqan::CharString> input_ids;
        std::vector<seqan::Dna5String> input_reads;
        std::vector<seqan::CharString> umi_ids;
        std::vector<seqan::Dna5String> umis;
        std::vector<seqan::Dna5String> compressed_umis;
        SparseGraphPtr umi_graph;
    };

    Input read_everything(const Params& params) {
        vector<seqan::CharString> input_ids;
        std::vector<seqan::Dna5String> input_reads;
        INFO("Reading records from " << params.reads_path);
        seqan::SeqFileIn reads_file(params.reads_path.c_str());
        readRecords(input_ids, input_reads, reads_file);
        INFO(input_ids.size() << " records read");

        vector<seqan::CharString> umi_ids;
        std::vector<seqan::Dna5String> umis;
        INFO("Reading UMI records from " << params.umi_uncompressed_path);
        seqan::SeqFileIn umi_file(params.umi_uncompressed_path.c_str());
        readRecords(umi_ids, umis, umi_file);
        INFO(umi_ids.size() << " UMIs read");
        VERIFY_MSG(umi_ids.size() == input_ids.size(), "The number of reads should be the same as the number of UMIs")

        std::vector<seqan::Dna5String> compressed_umis;
        INFO("Reading compressed UMIs from " << params.umi_compressed_path);
        std::vector<seqan::CharString> compressed_umi_ids;
        seqan::SeqFileIn umi_compressed_file(params.umi_compressed_path.c_str());
        readRecords(compressed_umi_ids, compressed_umis, umi_compressed_file);
        INFO(compressed_umis.size() << " compressed UMIs read");

        SparseGraphPtr umi_graph;
        INFO("Reading UMI graph from " << params.umi_graph_path);
        umi_graph = GraphReader(params.umi_graph_path).CreateGraph();
        INFO("Read graph with " << umi_graph->N() << " vertices and " << umi_graph->NZ() << " edges");

        return Input(input_ids, input_reads, umi_ids, umis, compressed_umis, umi_graph);
    }
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    const auto& input = read_everything(params);

    // needs uncompressed umis
    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(input.umis, umi_to_reads);

    // TODO: get rid of shared_ptrs to UMIs (bare UMIs even weigh less),
    // TODO: also share reads (but not individually),
    // TODO: use indices instead of map keys and values where possible
    std::unordered_map<Umi, UmiPtr> umi_ptr_by_umi;
    for (auto& umi_read : input.compressed_umis) {
        umi_ptr_by_umi[Umi(umi_read)] = std::make_shared<Umi>(umi_read);
    }

    std::vector<UmiPtr> compressed_umi_ptrs;
    for (auto& entry : umi_ptr_by_umi) {
        compressed_umi_ptrs.push_back(entry.second);
    }
    std::vector<Read> reads;
    for (size_t i = 0; i < input.input_reads.size(); i ++) {
        reads.emplace_back(input.input_reads[i], input.input_ids[i], i);
    }
    clusterer::ManyToManyCorrespondenceUmiToCluster<Read> initial_umis_to_clusters;
    for (auto& entry : umi_to_reads) {
        const auto& umi = umi_ptr_by_umi[entry.first];
        for (auto& read_idx : entry.second) {
            const auto& cluster = std::make_shared<clusterer::Cluster<Read>>(reads[read_idx], input.input_reads[read_idx], read_idx);
            initial_umis_to_clusters.add(umi, cluster);
        }
    }
    // This works 10-15 times slower than simple way (2 min vs 10 sec). Maybe because we don't try to glue to larger clusters first. Anyway, not a great problem at the moment.
    INFO("Clustering reads by hamming within single UMIs with threshold " << clusterer::ClusteringMode::hamming.threshold);
    const auto umi_to_clusters_hamm_inside_umi = clusterer::Clusterer<Read, clusterer::ReflexiveUmiPairsIterable>::cluster(
            clusterer::ClusteringMode::hamming, compressed_umi_ptrs, initial_umis_to_clusters,
            clusterer::ReflexiveUmiPairsIterable(compressed_umi_ptrs.size()));
    for (const auto& cluster : umi_to_clusters_hamm_inside_umi.toSet()) {
        VERIFY_MSG(umi_to_clusters_hamm_inside_umi.back(cluster).size() == 1, "We haven't united any reads across different UMIs yet.");
    }
    INFO(umi_to_clusters_hamm_inside_umi.toSize() << " clusters found");

    std::map<size_t, size_t> umis_of_size;
    std::map<size_t, size_t> clusters_in_umis_of_size;
    for (const auto& umi_ptr : compressed_umi_ptrs) {
        size_t size = umi_to_reads[*umi_ptr].size();
        umis_of_size[size] ++;
        clusters_in_umis_of_size[size] += umi_to_clusters_hamm_inside_umi.forth(umi_ptr).size();
    }
    for (const auto& entry : umis_of_size) {
        size_t size = entry.first;
        INFO("Size: " << size << ", UMIs: " << entry.second << ", clusters: " << clusters_in_umis_of_size[size] << ", compression: "
                     << 100.0 * static_cast<double>(clusters_in_umis_of_size[size]) / static_cast<double>(size * entry.second) << "%");
    }

    INFO("Uniting read clusters for adjacent UMIs");
    const auto umi_to_clusters_hamm_adj_umi = clusterer::Clusterer<Read, clusterer::GraphUmiPairsIterable>::cluster(
            clusterer::ClusteringMode::hamming, compressed_umi_ptrs, /*initial_umis_to_clusters*/umi_to_clusters_hamm_inside_umi,
            clusterer::GraphUmiPairsIterable(input.umi_graph));
    INFO(umi_to_clusters_hamm_adj_umi.toSize() << " clusters found");

    INFO("Clustering reads by edit distance within single UMIs with threshold " << clusterer::ClusteringMode::edit.threshold);
    const auto umi_to_clusters_edit_inside_umi = clusterer::Clusterer<Read, clusterer::ReflexiveUmiPairsIterable>::cluster(
            clusterer::ClusteringMode::edit, compressed_umi_ptrs, umi_to_clusters_hamm_adj_umi,
            clusterer::ReflexiveUmiPairsIterable(compressed_umi_ptrs.size()));
    INFO(umi_to_clusters_edit_inside_umi.toSize() << " clusters found");

    INFO("Uniting read clusters for adjacent UMIs");
    const auto umi_to_clusters_edit_adj_umi = clusterer::Clusterer<Read, clusterer::GraphUmiPairsIterable>::cluster(
            clusterer::ClusteringMode::edit, compressed_umi_ptrs, umi_to_clusters_edit_inside_umi,
            clusterer::GraphUmiPairsIterable(input.umi_graph));
    INFO(umi_to_clusters_edit_adj_umi.toSize() << " clusters found");

    INFO("Uniting clusters with identical centers");
    const auto umi_to_clusters_same_centers = clusterer::Clusterer<Read, clusterer::GraphUmiPairsIterable>::uniteByCenters(
            umi_to_clusters_edit_adj_umi);
    INFO(umi_to_clusters_same_centers.toSize() << " clusters found");

//    INFO("Uniting ignoring UMIs");
//    const auto& umi_to_clusters_global = clusterer::Clusterer<Read, clusterer::FullGraphUmiPairsIterable>::cluster(
//            clusterer::ClusteringMode::hamming, compressed_umi_ptrs, umi_to_clusters_edit_adj_umi,
//            clusterer::FullGraphUmiPairsIterable(compressed_umi_ptrs.size()));
//    INFO(umi_to_clusters_global.toSize() << " clusters found");

    INFO("Saving intermediate repertoire to output directory " << params.output_dir);
    clusterer::write_clusters_and_correspondence<Read>(umi_to_clusters_same_centers, reads, params.output_dir, params.save_clusters);
    INFO("Saving finished");
    // unite close reads with different UMIs: graph is needed anyway; then either metis clustering, or continue custom techniques

    return 0;
}
