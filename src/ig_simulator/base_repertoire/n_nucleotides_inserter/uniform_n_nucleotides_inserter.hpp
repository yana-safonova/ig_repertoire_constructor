//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "abstract_n_nucleotides_inserter.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class UniformNNucleotidesInserter final : public AbstractNNucleotidesInserter {
private:
    const size_t max_vj_insertion;
    const size_t max_vd_insertion;
    const size_t max_dj_insertion;

    seqan::Dna5String RandDna5Str(size_t size) const;

public:
    explicit UniformNNucleotidesInserter(
        const NNucleotidesInserterParams::UniformInserterParams config):
                max_vj_insertion(config.max_vj_insertion),
                max_vd_insertion(config.max_vd_insertion),
                max_dj_insertion(config.max_dj_insertion)
    { }

    seqan::Dna5String GetVJInsertion() const override;
    seqan::Dna5String GetVDInsertion() const override;
    seqan::Dna5String GetDJInsertion() const override;
};

} // End namespace ig_simulator
