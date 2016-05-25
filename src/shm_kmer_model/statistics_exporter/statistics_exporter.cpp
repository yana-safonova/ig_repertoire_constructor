//
// Created by Andrew Bzikadze on 5/24/16.
//

#include <fstream>

#include "seqan/basic.h"
#include "statistics_exporter.hpp"

void StatisticsExporter::export_statistics(const MutationsStatistics &statistics) const {
    std::ofstream out(output_filename_);

    out << "k-mer";
    auto alphabet_size = seqan::ValueSize<seqan::Dna>::VALUE;
    for (size_t i = 0; i < alphabet_size; ++i)
        out << separator << seqan::Dna(i);
    out << "\n";

    std::vector<std::string> kmers;
    for (auto it = statistics.cbegin(); it != statistics.cend(); ++it)
        kmers.push_back(it -> first);
    std::sort(kmers.begin(), kmers.end());

    for (auto it = kmers.begin(); it != kmers.end(); ++it) {
        auto& vec = statistics.at(*it);
        out << *it;
        for (auto it2=vec.begin(); it2 != vec.end(); ++it2)
            out <<separator << *it2 ;
        out << "\n";
    }
}
