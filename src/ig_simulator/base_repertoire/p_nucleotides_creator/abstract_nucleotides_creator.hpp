//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include <cstdlib>
#include <memory>

namespace ig_simulator {

class AbstractPNucleotidesCreator {
public:
    virtual size_t CreateInVGene() const = 0;
    virtual size_t CreateInDGeneLeft() const = 0;
    virtual size_t CreateInDGeneRight() const = 0;
    virtual size_t CreateInJGene() const = 0;

    AbstractPNucleotidesCreator() = default;
    AbstractPNucleotidesCreator(const AbstractPNucleotidesCreator&) = delete;
    AbstractPNucleotidesCreator(AbstractPNucleotidesCreator&&) = delete;
    AbstractPNucleotidesCreator& operator=(const AbstractPNucleotidesCreator&) = delete;
    AbstractPNucleotidesCreator& operator=(AbstractPNucleotidesCreator&&) = delete;

    virtual ~AbstractPNucleotidesCreator() { }
};

using AbstractPNucleotidesCreatorCPtr = std::unique_ptr<const AbstractPNucleotidesCreator>;

} // End namespace ig_simulator
