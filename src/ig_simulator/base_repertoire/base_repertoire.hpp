//
// Created by Andrew Bzikadze on 4/2/17.
//

#pragma once

#include "metaroot_cluster/metaroot_cluster.hpp"

namespace ig_simulator {

using BaseRepertoire = std::vector<MetarootCluster>;

std::ostream& operator<<(std::ostream&, const BaseRepertoire&);


} // End namespace ig_simulator
