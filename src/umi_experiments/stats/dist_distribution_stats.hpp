#pragma once

#include "umi_utils.hpp"

class DistDistributionStats {
public:
    DistDistributionStats(std::map<size_t, std::map<size_t, size_t>> hamming_dist_distribution, std::map<size_t, std::map<size_t, size_t>> min_hamming_dist_distribution,
                          std::map<size_t, std::map<size_t, size_t>> sw_dist_distribution, std::map<size_t, std::map<size_t, size_t>> min_sw_dist_distribution)
            : hamming_dist_distribution_(hamming_dist_distribution), min_hamming_dist_distribution_(min_hamming_dist_distribution),
              sw_dist_distribution_(sw_dist_distribution), min_sw_dist_distribution_(min_sw_dist_distribution) {};

    std::vector<size_t> GetSizes();

    std::string ToString(size_t umi_size);

    static DistDistributionStats GetStats(const std::vector<seqan::Dna5String>& input_reads, const std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads, unsigned thread_count);

private:
    std::map<size_t, std::map<size_t, size_t>> hamming_dist_distribution_;
    std::map<size_t, std::map<size_t, size_t>> min_hamming_dist_distribution_;
    std::map<size_t, std::map<size_t, size_t>> sw_dist_distribution_;
    std::map<size_t, std::map<size_t, size_t>> min_sw_dist_distribution_;
};
