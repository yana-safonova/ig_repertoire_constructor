//
// Created by Andrew Bzikadze on 4/2/17.
//

#include "base_repertoire.hpp"

namespace ig_simulator {

void print_base_repertoire(const BaseRepertoire& base_repertoire, std::ostream& fasta, std::ostream& info) {
    size_t id = 0;
    for (const auto& cluster : base_repertoire) {
        fasta << ">forest_" << id << "_multiplicity_" << cluster.Multiplicity() << '\n';
        fasta << cluster.MetarootPtr()->Sequence() << '\n';

        info << "Index (zero-based): " << id << " / " << base_repertoire.size() - 1
             << " (" << base_repertoire.size() << ") "
             << '\n' << *(cluster.MetarootPtr())
             <<"***************************************************************************\n\n";

        id++;
    }
}

} // End namespapce ig_simulator
