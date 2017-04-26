//
// Created by Andrew Bzikadze on 3/29/17.
//

#include "base_repertoire_simulator.hpp"

namespace ig_simulator {

BaseRepertoire BaseRepertoireSimulator::Simulate(size_t size) {
    BaseRepertoire repertoire;
    repertoire.reserve(size);

    size_t productive_size = static_cast<size_t> (static_cast<double>(size) * productive_part);
    size_t i = 0;
    while(i < productive_size) {
        MetarootCluster cluster{metaroot_creator_p->Createroot(),
                                multiplicity_creator_p->RandomMultiplicity()};
        if (cluster.MetarootPtr()->IsProductive()) {
            repertoire.emplace_back(std::move(cluster));
            i++;
        }
    }

    for(; i < size; ++i) {
        repertoire.emplace_back(metaroot_creator_p->Createroot(),
                                multiplicity_creator_p->RandomMultiplicity());
    }
    return repertoire;
}

} // End namespace ig_simulator
