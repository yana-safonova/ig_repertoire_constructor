#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <seqan/seq_io.h>
#include "utils.hpp"
#include "umi_utils.hpp"
#include "../graph_utils/graph_io.hpp"
#include "clusterer.hpp"

namespace {
    struct Params {
        std::string reads_path;
        std::string umi_uncompressed_path;
        std::string umi_compressed_path;
        std::string umi_graph_path;
        std::string output_dir;
        bool detect_chimeras;
        bool save_clusters;
        size_t num_threads;
        size_t clustering_threshold;
        bool output_intermediate;
    };

    bool read_args(int argc, const char* const* argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("reads,r", po::value<std::string>(&params.reads_path)->required(), "input file with reads")
                ("umi-uncompressed,u", po::value<std::string>(&params.umi_uncompressed_path)->required(), "file with UMI records extracted (not compressed)")
                ("umi-compressed,c", po::value<std::string>(&params.umi_compressed_path)->required(), "file with UMI records extracted (compressed)")
                ("graph,g", po::value<std::string>(&params.umi_graph_path)->required(), "file with UMI graph")
                ("output,o", po::value<std::string>(&params.output_dir)->default_value(""), "output directory path")
                ("detect-chimeras,k", po::value<bool>(&params.detect_chimeras)->default_value(false), "detect chimeras after clustering, may take significant amount of time")
                ("save-clusters,s", po::value<bool>(&params.save_clusters)->default_value(false), "save clusters by UMI")
                ("threads,t", po::value<size_t >(&params.num_threads)->default_value(1), "number of threads to use")
                ("clustering-thr,d", po::value<size_t >(&params.clustering_threshold)->default_value(20), "threshold distance to unite clusters")
                ("debug-stages,b", po::value<bool>(&params.output_intermediate)->default_value(false), "output repertoire after each step")
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

int main(int argc, const char* const* argv) {
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
    std::vector<UmiPtr> compressed_umi_ptrs;
    for (auto& umi_read : input.compressed_umis) {
        const auto& umi_ptr = std::make_shared<Umi>(umi_read);
        compressed_umi_ptrs.push_back(umi_ptr);
        umi_ptr_by_umi[Umi(umi_read)] = umi_ptr;
    }

    std::vector<Read> reads;
    for (size_t i = 0; i < input.input_reads.size(); i ++) {
        reads.emplace_back(input.input_reads[i], input.input_ids[i], i);
    }
    clusterer::Clusterer<Read> clusterer(umi_to_reads, umi_ptr_by_umi, reads);
    // This works 10-15 times slower than simple way (2 min vs 10 sec). Maybe because we don't try to glue to larger clusters first. Anyway, not a great problem at the moment.
    const clusterer::ReadDist& hamming_dist = clusterer::ClusteringMode::bounded_hamming_dist(params.clustering_threshold);
    const auto hamming_dist_checker = clusterer::ClusteringMode::clusters_close_by_min(hamming_dist, params.clustering_threshold);
    INFO("Clustering reads by hamming within single UMIs with threshold " << params.clustering_threshold);
    clusterer.cluster(hamming_dist_checker, compressed_umi_ptrs, clusterer::ReflexiveUmiPairsIterable(compressed_umi_ptrs.size()));
    for (const auto& cluster : clusterer.getCurrentUmiToCluster().toSet()) {
        VERIFY_MSG(clusterer.getCurrentUmiToCluster().back(cluster).size() == 1, "We haven't united any reads across different UMIs yet.");
    }
    INFO(clusterer.getCurrentUmiToCluster().toSize() << " clusters found");

    clusterer.write_clusters_and_correspondence(params.output_dir, "_1_hamming", params.save_clusters, params.output_intermediate);


//    clusterer::report_non_major_umi_groups(umi_to_clusters_hamm_inside_umi, params.output_dir + "/non_major.csv");

//    std::map<size_t, size_t> umis_of_size;
//    std::map<size_t, size_t> clusters_in_umis_of_size;
//    for (const auto& umi_ptr : compressed_umi_ptrs) {
//        size_t size = umi_to_reads[*umi_ptr].size();
//        umis_of_size[size] ++;
//        clusters_in_umis_of_size[size] += umi_to_clusters_hamm_inside_umi.forth(umi_ptr).size();
//    }
//    for (const auto& entry : umis_of_size) {
//        size_t size = entry.first;
//        INFO("Size: " << size << ", UMIs: " << entry.second << ", clusters: " << clusters_in_umis_of_size[size] << ", compression: "
//                     << 100.0 * static_cast<double>(clusters_in_umis_of_size[size]) / static_cast<double>(size * entry.second) << "%");
//    }

    clusterer.print_umi_split_stats();

    INFO("Uniting read clusters for adjacent UMIs");
    clusterer.cluster(hamming_dist_checker, compressed_umi_ptrs, clusterer::GraphUmiPairsIterable(input.umi_graph));
    INFO(clusterer.getCurrentUmiToCluster().toSize() << " clusters found");
//    size_t hamm_corrected_reads = clusterer::count_reads_with_corrected_umi(umi_to_clusters_hamm_inside_umi, umi_to_clusters_hamm_adj_umi);
//    INFO(hamm_corrected_reads << " reads have UMI corrected for hamming dist.");

    clusterer.write_clusters_and_correspondence(params.output_dir, "_2_hamming_ngh", params.save_clusters, params.output_intermediate);


    const clusterer::ReadDist& edit_dist = clusterer::ClusteringMode::bounded_edit_dist(params.clustering_threshold, params.clustering_threshold);
    const auto edit_dist_checker = clusterer::ClusteringMode::clusters_close_by_min(edit_dist, params.clustering_threshold);
    INFO("Clustering reads by edit distance within single UMIs with threshold " << params.clustering_threshold);
    clusterer.cluster(edit_dist_checker, compressed_umi_ptrs, clusterer::ReflexiveUmiPairsIterable(compressed_umi_ptrs.size()));
    INFO(clusterer.getCurrentUmiToCluster().toSize() << " clusters found");

    clusterer.write_clusters_and_correspondence(params.output_dir, "_3_edit", params.save_clusters, params.output_intermediate);


    INFO("Uniting read clusters for adjacent UMIs");
    clusterer.cluster(edit_dist_checker, compressed_umi_ptrs, clusterer::GraphUmiPairsIterable(input.umi_graph));
    INFO(clusterer.getCurrentUmiToCluster().toSize() << " clusters found");

    clusterer.write_clusters_and_correspondence(params.output_dir, "_4_edit_ngh", params.save_clusters, params.output_intermediate);


    clusterer.report_length_differences(params.output_dir);

    if (params.detect_chimeras) {
        clusterer.report_non_major_umi_groups_sw(params.output_dir + "/non_major.csv",
                                                 params.output_dir + "/left_graph.graph",
                                                 params.output_dir + "/right_graph.graph",
                                                 params.output_dir + "/chimeras.txt",
                                                 params.output_dir + "/umi_chimeras.txt",
                                                 params.num_threads);

        clusterer.write_clusters_and_correspondence(params.output_dir, "_5_chimeras", params.save_clusters, params.output_intermediate);
    }
    clusterer.write_clusters_and_correspondence(params.output_dir, "", params.save_clusters, true, "", "intermediate_repertoire_close_umi.rcm");

//    size_t edit_corrected_reads = clusterer::count_reads_with_corrected_umi(umi_to_clusters_edit_inside_umi, umi_to_clusters_edit_adj_umi);
//    INFO(edit_corrected_reads << " reads have UMI corrected for edit dist.");
//    INFO(hamm_corrected_reads + edit_corrected_reads << " reads total have UMI corrected.");
    size_t reads_with_corrected_umis = clusterer::count_reads_with_corrected_umi(clusterer.getCurrentUmiToCluster());
    INFO(reads_with_corrected_umis << " reads total have UMI corrected.");

    INFO("Uniting clusters with identical centers");
    clusterer.uniteByCenters();
    INFO(clusterer.getCurrentUmiToCluster().toSize() << " clusters found");

//    INFO("Uniting ignoring UMIs");
//    const auto& umi_to_clusters_global = clusterer::Clusterer<Read, clusterer::FullGraphUmiPairsIterable>::cluster(
//            clusterer::ClusteringMode::hamming, compressed_umi_ptrs, umi_to_clusters_edit_adj_umi,
//            clusterer::FullGraphUmiPairsIterable(compressed_umi_ptrs.size()));
//    INFO(umi_to_clusters_global.toSize() << " clusters found");

    INFO("Saving intermediate repertoire to output directory " << params.output_dir);
    if (params.output_intermediate) {
        clusterer.write_clusters_and_correspondence(params.output_dir, "_6_unite_equal", params.save_clusters, params.output_intermediate);
    } else {
        clusterer.write_clusters_and_correspondence(params.output_dir, "", params.save_clusters, true);
    }
    INFO("Saving finished");

    return 0;
}
