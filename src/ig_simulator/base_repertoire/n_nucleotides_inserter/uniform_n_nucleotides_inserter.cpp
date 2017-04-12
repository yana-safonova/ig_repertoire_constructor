//
// Created by Andrew Bzikadze on 3/20/17.
//

#include <seqan/file.h>
#include "uniform_n_nucleotides_inserter.hpp"
#include "simulation_routines.hpp"

using seqan::Dna5String;

namespace ig_simulator {

Dna5String UniformNNucleotidesInserter::RandDna5Str(size_t size) const {
    auto RandomNucleotide = []() -> char {
      return "ACGT"[random_index(0, 3)];
    };

    std::vector<char> v_str(size);
    for (auto & nucl : v_str) {
        nucl = RandomNucleotide();
    }
    return Dna5String(std::string(v_str.begin(), v_str.end()));
}

Dna5String UniformNNucleotidesInserter::GetVJInsertion() const { return RandDna5Str(random_index(0, max_vj_insertion)); }
Dna5String UniformNNucleotidesInserter::GetVDInsertion() const { return RandDna5Str(random_index(0, max_vd_insertion)); }
Dna5String UniformNNucleotidesInserter::GetDJInsertion() const { return RandDna5Str(random_index(0, max_dj_insertion)); }

} // End namespace ig_simulator
