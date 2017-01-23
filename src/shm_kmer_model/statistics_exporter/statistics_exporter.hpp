//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include <string>

#include "../statistics_estimator/mutation_statistics.hpp"

namespace shm_kmer_matrix_estimator {

class StatisticsExporter {
private:
    const char separator = ';';

public:
    void export_statistics(const std::string &output_filename, const MutationsStatistics &) const;
};

} // End namespace shm_kmer_matrix_estimator
