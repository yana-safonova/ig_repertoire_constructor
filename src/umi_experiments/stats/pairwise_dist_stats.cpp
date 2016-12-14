#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <utils/io.hpp>
#include <clusterer.hpp>

void calc_and_write(int threads, const std::vector<seqan::Dna5String>& reads, size_t max_dist, std::ofstream& ofs) {
    omp_set_num_threads(threads);

    const size_t size = reads.size();
    const auto& get_dist = clusterer::ClusteringMode::bounded_edit_dist(max_dist, max_dist, false);
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t i = 0; i < size; i ++) {
        std::vector<size_t> dist(size - i - 1);

        for (size_t j = i + 1; j < size; j ++) {
            dist[j - i - 1] = get_dist(reads[i], reads[j]);
        }
        SEQAN_OMP_PRAGMA(critical)
        {
            std::for_each(dist.begin(), dist.end(), [&ofs, &max_dist](size_t d) { if (d <= max_dist) ofs << d << "\n"; });
        }
        std::cout << i << "\n" << clock() << std::endl;
    }
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;

    const std::string repertoire_file(argv[1]);
    std::vector<seqan::CharString> ids;
    std::vector<seqan::Dna5String> reads;
    read_seqan_records(repertoire_file, ids, reads);

    const std::string dist_output_file(argv[2]);
    std::ofstream ofs(dist_output_file);
    const int threads = std::stoi(argv[3]);
    const size_t max_dist = std::stoull(argv[4]);
    calc_and_write(threads, reads, max_dist, ofs);
}
