#include <unordered_map>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <utils/io.hpp>
#include <verify.hpp>
#include <clusterer.hpp>
#include "../../ig_tools/utils/string_tools.hpp"

std::unordered_map<std::string, seqan::Dna5String> get_id_to_read_map(const std::string reads_file) {
    std::vector<seqan::CharString> ids;
    std::vector<seqan::Dna5String> reads;
    read_seqan_records(reads_file, ids, reads);
    std::unordered_map<std::string, seqan::Dna5String> id2read;
    for (size_t i = 0; i < ids.size(); i ++) {
        id2read[seqan_string_to_string(ids[i])] = reads[i];
    }
    return id2read;
}

std::vector<std::vector<std::string>> get_clusters_from_rcm(const std::string& rcm_file) {
    INFO("Reading RCM from " << rcm_file);
    std::unordered_map<std::string, std::vector<std::string>> clusters;
    std::ifstream rcm(rcm_file);
    std::string line;
    while (std::getline(rcm, line)) {
        std::vector<std::string> parts = split(line, '\t');
        VERIFY(parts.size() <= 2 && !parts.empty());
        if (parts.size() == 1) continue;
        clusters[parts[1]].push_back(parts[0]);
    }
    std::vector<std::vector<std::string>> result(clusters.size());
    size_t cur = 0;
    for (const auto& entry : clusters) {
        result[cur ++] = entry.second;
    }
    INFO(result.size() << " records read");
    return result;
}

int main(int, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();

    const std::string reads_file(argv[1]);
    const std::string rcm_file(argv[2]);
    const std::string dist_output_file(argv[3]);
    const int threads = std::stoi(argv[4]);
    const size_t max_dist = std::stoull(argv[5]);
    
    auto id2read = get_id_to_read_map(reads_file);
    const auto clusters = get_clusters_from_rcm(rcm_file);
    
    const auto get_dist = clusterer::ClusteringMode::bounded_edit_dist(max_dist, max_dist, false);
    std::ofstream ofs(dist_output_file);
    omp_set_num_threads(threads);
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t cluster = 0; cluster < clusters.size(); cluster ++) {
        INFO(cluster);
        const size_t size = clusters[cluster].size();
        size_t cur = 0;
        std::vector<size_t> dists(size * (size - 1) / 2);
        for (size_t i = 0; i < size; i ++) {
            for (size_t j = 0; j < i; j ++) {
                dists[cur ++] = get_dist(id2read[clusters[cluster][i]], id2read[clusters[cluster][j]]);
            }
        }

        SEQAN_OMP_PRAGMA(critical)
        if (size > 1) {
            ofs << size << "\n";
            for (size_t dist : dists) {
                ofs << dist << "\n";
            }
        }
    }
}
