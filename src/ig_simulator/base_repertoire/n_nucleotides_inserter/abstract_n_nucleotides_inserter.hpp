//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include <cstdlib>
#include <memory>
#include <seqan/seq_io.h>

namespace ig_simulator {

class AbstractNNucleotidesInserter {
public:
    virtual seqan::Dna5String GetVJInsertion() const = 0;
    virtual seqan::Dna5String GetVDInsertion() const = 0;
    virtual seqan::Dna5String GetDJInsertion() const = 0;

    AbstractNNucleotidesInserter() = default;
    AbstractNNucleotidesInserter(const AbstractNNucleotidesInserter&) = delete;
    AbstractNNucleotidesInserter(AbstractNNucleotidesInserter&&) = delete;
    AbstractNNucleotidesInserter& operator=(const AbstractNNucleotidesInserter&) = delete;
    AbstractNNucleotidesInserter& operator=(AbstractNNucleotidesInserter&&) = delete;

    virtual ~AbstractNNucleotidesInserter() { }
};

using AbstractNNucleotidesInserterCPtr = std::unique_ptr<const AbstractNNucleotidesInserter>;

} // End namespace ig_simulator
