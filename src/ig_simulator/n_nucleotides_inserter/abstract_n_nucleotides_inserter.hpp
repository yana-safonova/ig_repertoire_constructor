//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include <cstdlib>
#include <memory>

namespace ig_simulator {

class AbstractNNucleotidesInserter {
public:
    virtual size_t GetVJInsertion() const = 0;
    virtual size_t GetVDInsertion() const = 0;
    virtual size_t GetDJInsertion() const = 0;

    virtual ~AbstractNNucleotidesInserter() { }
};

using AbstractNNucleotidesInserterPtr = std::unique_ptr<AbstractNNucleotidesInserter>;
} // End namespace ig_simulator
