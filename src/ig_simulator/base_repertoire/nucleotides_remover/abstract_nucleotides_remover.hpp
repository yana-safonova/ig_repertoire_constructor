//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <cstdlib>
#include <memory>

namespace ig_simulator {

class AbstractNucleotidesRemover {
public:
    virtual size_t RemoveInVGene()      const = 0;
    virtual size_t RemoveInDGeneLeft()  const = 0;
    virtual size_t RemoveInDGeneRight() const = 0;
    virtual size_t RemoveInJGene()      const = 0;

    AbstractNucleotidesRemover() = default;
    AbstractNucleotidesRemover(const AbstractNucleotidesRemover&) = delete;
    AbstractNucleotidesRemover(AbstractNucleotidesRemover&&) = delete;
    AbstractNucleotidesRemover& operator=(const AbstractNucleotidesRemover&) = delete;
    AbstractNucleotidesRemover& operator=(AbstractNucleotidesRemover&&) = delete;

    virtual ~AbstractNucleotidesRemover() { }
};

using AbstractNucleotidesRemoverCPtr = std::unique_ptr<const AbstractNucleotidesRemover>;

} // End namespace ig_simulator
