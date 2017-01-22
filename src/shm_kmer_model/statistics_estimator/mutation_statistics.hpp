//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <seqan/file.h>

class MutationsStatistics {
private:
    using statistics_type = std::unordered_map<std::string, std::vector<unsigned int>>;
    statistics_type statistics;

public:
    explicit MutationsStatistics(unsigned int kmer_len = unsigned());

    statistics_type::iterator begin() noexcept { return statistics.begin(); }
    statistics_type::const_iterator begin() const noexcept { return statistics.begin(); }
    statistics_type::const_iterator cbegin() const noexcept { return statistics.cbegin(); }

    statistics_type::iterator end() noexcept { return statistics.end(); }
    statistics_type::const_iterator end() const noexcept { return statistics.end(); }
    statistics_type::const_iterator cend() const noexcept { return statistics.cend(); }

    statistics_type::iterator find(const std::string &str) { return statistics.find(str); }
    statistics_type::const_iterator find(const std::string &str) const { return statistics.find(str); }

    std::vector<unsigned int> &operator[](const std::string &str) { return statistics[str]; }
    std::vector<unsigned int> &operator[](std::string &str) { return statistics[str]; }

    std::vector<unsigned int> &at(const std::string &str) { return statistics.at(str); }
    const std::vector<unsigned int> &at(const std::string &str) const { return statistics.at(str); }

private:
    void generate_kmer_keys(std::vector<std::string> &kmers, const unsigned int kmer_len) const;

    void generate_kmer_keys_(std::vector<std::string> &kmers,
                             const unsigned int kmer_len,
                             unsigned int curr_len,
                             std::string kmer) const;
};