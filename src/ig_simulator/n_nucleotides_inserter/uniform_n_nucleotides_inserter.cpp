//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "uniform_n_nucleotides_inserter.hpp"
#include "random_index.hpp"

namespace ig_simulator {

size_t UniformNNucleotidesInserter::GetVJInsertion() const { return random_index(0, max_vj_insertion); }
size_t UniformNNucleotidesInserter::GetVDInsertion() const { return random_index(0, max_vd_insertion); }
size_t UniformNNucleotidesInserter::GetDJInsertion() const { return random_index(0, max_dj_insertion); }

} // End namespace ig_simulator
