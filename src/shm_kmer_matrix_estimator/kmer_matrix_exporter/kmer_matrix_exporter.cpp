//
// Created by Andrew Bzikadze on 5/24/16.
//

#include <fstream>

#include "seqan/basic.h"
#include "kmer_matrix_exporter.hpp"
#include "kmer_utils/kmer_utils.hpp"

namespace shm_kmer_matrix_estimator {

void KmerMatrixExporter::export_statistics(const std::string &output_filename,
                                           const KmerMatrix &statistics) const {
    std::ofstream out(output_filename);

    auto alphabet_size = seqan::ValueSize<seqan::Dna>::VALUE;
    for (size_t i = 0; i < alphabet_size; ++i)
        out << separator << seqan::Dna(i);
    out << "\n";

    // TODO remove magic const
    std::vector<std::string> kmers(KmerUtils::GenerateAllKmersFixedLength(5));

    for (const auto& kmer: kmers) {
        const std::array<unsigned int, 4>& vec(statistics.at(kmer));
        out << kmer;
        for (const auto& freq: vec)
            out << separator << freq;
        out << "\n";
    }
}

} // End namespace shm_kmer_matrix_estimator
