//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include <string>

#include "../statistics_estimator/mutation_statistics.hpp"

class StatisticsExporter {
private:
    std::string output_filename_;
    const char separator = ';';

public:
    explicit StatisticsExporter(const std::string &output_filename) :
        output_filename_(output_filename) { }

    void export_statistics(const MutationsStatistics &) const;
};
