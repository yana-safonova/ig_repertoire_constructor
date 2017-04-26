//
// Created by Andrew Bzikadze on 4/2/17.
//

#pragma once

#include "metaroot_cluster/metaroot_cluster.hpp"

namespace ig_simulator {

using BaseRepertoire = std::vector<MetarootCluster>;

void print_base_repertoire(const BaseRepertoire& base_repertoire, std::ostream& fasta, std::ostream& info);

} // End namespace ig_simulator
