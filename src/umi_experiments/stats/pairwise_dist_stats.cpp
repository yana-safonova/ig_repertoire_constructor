#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <utils/io.hpp>
#include <clusterer.hpp>

void calc_and_write(int threads, const std::vector<seqan::Dna5String>& reads, size_t max_dist, std::ofstream& ofs) {
    omp_set_num_threads(threads);

    const size_t size = reads.size();
    std::vector<size_t> dists(max_dist + 2);
    const auto& get_dist = clusterer::ClusteringMode::bounded_edit_dist(max_dist, max_dist, false);
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t i = 0; i < size; i ++) {
        std::vector<size_t> local_dists(max_dist + 2);

        for (size_t j = i + 1; j < size; j ++) {
            local_dists[get_dist(reads[i], reads[j])] ++;
        }
        SEQAN_OMP_PRAGMA(critical)
        for (size_t d = 0; d <= max_dist + 1; d ++) {
            dists[d] += local_dists[d];
        }
        std::cout << i << "\n" << clock() << std::endl;
    }
    std::for_each(dists.begin(), dists.end(), [&ofs](size_t d) { ofs << d << "\n"; });
}

int main(int, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();

    const std::string reads_file(argv[1]);
    std::vector<seqan::CharString> ids;
    std::vector<seqan::Dna5String> reads;
    read_seqan_records(reads_file, ids, reads);

    const std::string dist_output_file(argv[2]);
    std::ofstream ofs(dist_output_file);
    const int threads = std::stoi(argv[3]);
    const size_t max_dist = std::stoull(argv[4]);
    calc_and_write(threads, reads, max_dist, ofs);
}
