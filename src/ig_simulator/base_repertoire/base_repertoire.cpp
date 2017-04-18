//
// Created by Andrew Bzikadze on 4/2/17.
//

#include "base_repertoire.hpp"

namespace ig_simulator {

std::ostream &operator<<(std::ostream &out, const BaseRepertoire &base_repertoire) {
    size_t id = 0;
    for (const auto& cluster : base_repertoire) {
        out << ">forest_" << id++ << "_multiplicity_" << cluster.Multiplicity() << '\n';
        out << cluster.MetarootPtr()->Sequence() << '\n';
    }
    return out;
}

} // End namespapce ig_simulator
