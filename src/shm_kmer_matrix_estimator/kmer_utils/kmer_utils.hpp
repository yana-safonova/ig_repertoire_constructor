//
// Created by Andrew Bzikadze on 11/9/16.
//

#pragma once

#include <cstdlib>
#include <string>
#include <vector>

namespace shm_kmer_matrix_estimator {

class KmerUtils {
private:
    static void GenerateAllKmersFixedLength(const size_t kmer_len,
                                            const size_t position,
                                            std::string& current_kmer,
                                            std::vector<std::string>& kmer_vector);

public:
    static size_t GetIndexNucl(const char nucl);
    static size_t GetIndexByKmer(const std::string& kmer);
    static std::vector<std::string> GenerateAllKmersFixedLength(const size_t kmer_len);
    static size_t GetMutationIndexByKmerAndNucl(const std::string& kmer, const char nucl);
};

} // End namespace shm_kmer_matrix_estimator