//
// Created by Andrew Bzikadze on 3/29/17.
//

#include "base_repertoire_simulator.hpp"

namespace ig_simulator {

BaseRepertoire BaseRepertoireSimulator::Simulate(size_t size) {
    BaseRepertoire repertoire;
    repertoire.reserve(size);

    for (size_t i = 0; i < size; ++i) {
        repertoire.emplace_back(metaroot_creator_p->Createroot(),
                                multiplicity_creator_p->RandomMultiplicity());
    }

    return repertoire;
}

} // End namespace ig_simulator
