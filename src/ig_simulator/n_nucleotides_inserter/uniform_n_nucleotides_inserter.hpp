//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "abstract_n_nucleotides_inserter.hpp"

namespace ig_simulator {

class UniformNNucleotidesInserter : public AbstractNNucleotidesInserter {
private:
    const size_t max_vj_insertion = 10;
    const size_t max_vd_insertion = 10;
    const size_t max_dj_insertion = 10;

public:
    virtual size_t GetVJInsertion() const override;
    virtual size_t GetVDInsertion() const override;
    virtual size_t GetDJInsertion() const override;

    UniformNNucleotidesInserter() { }
};

} // End namespace ig_simulator
