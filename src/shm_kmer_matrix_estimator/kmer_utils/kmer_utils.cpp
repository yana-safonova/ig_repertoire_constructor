//
// Created by Andrew Bzikadze on 11/9/16.
//

#include "seqan/basic.h"
#include "verify.hpp"

#include "kmer_utils.hpp"

namespace shm_kmer_matrix_estimator {

size_t KmerUtils::GetIndexNucl(const char nucl) {
    return static_cast<size_t>(seqan::ordValue(static_cast<seqan::Dna>(nucl)));
}

size_t KmerUtils::GetIndexByKmer(const std::string &kmer) {

    size_t index_kmer = 0;

    size_t i = kmer.length();
    size_t pow = 1;
    while (i--) {
        size_t index_nucl = GetIndexNucl(kmer[i]);
        index_kmer += pow * index_nucl;
        pow *= seqan::ValueSize<seqan::Dna>::VALUE;
    }
    return index_kmer;
}

void KmerUtils::GenerateAllKmersFixedLength(const size_t kmer_len,
                                            const size_t position,
                                            std::string &current_kmer,
                                            std::vector<std::string> &kmer_vector) {
    if (position == kmer_len) {
        kmer_vector.push_back(current_kmer);
        return;
    }
    auto alphabet_size = seqan::ValueSize<seqan::Dna>::VALUE;
    for (size_t i = 0; i < alphabet_size; ++i) {
        current_kmer[position] = seqan::Dna(i);
        GenerateAllKmersFixedLength(kmer_len, position + 1, current_kmer, kmer_vector);
    }
}

std::vector<std::string> KmerUtils::GenerateAllKmersFixedLength(const size_t kmer_len) {
    std::vector<std::string> kmers;
    std::string current_kmer(kmer_len, ' ');
    GenerateAllKmersFixedLength(kmer_len, 0, current_kmer, kmers);
    return kmers;
}


size_t KmerUtils::GetMutationIndexByKmerAndNucl(const std::string &kmer, const char nucl) {
    const char center_kmer_nucl = kmer[kmer.size() / 2];
    VERIFY_MSG(center_kmer_nucl != nucl, "Mutation of central nucl to the same nucl");
    size_t index_nucl(GetIndexNucl(nucl));
    size_t index_center_kmer_nucl(GetIndexNucl(center_kmer_nucl));
    return index_nucl - (index_nucl > index_center_kmer_nucl);
}

} // End namespace shm_kmer_matrix_estimator
