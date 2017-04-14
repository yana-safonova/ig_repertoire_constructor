//
// Created by Andrew Bzikadze on 4/14/17.
//

#pragma once

#include <fstream>
#include <clonal_trees/tree/tree.hpp>
#include <clonal_trees/forest/forest.hpp>

namespace ig_simulator {

void TreeExport(const Tree& tree, size_t forest_ind, size_t tree_ind, std::ostream& full, std::ostream& included);
void ForestExporter(const Forest& forest, size_t forest_ind, std::ostream& full, std::ostream& included);
void ForestStorageExporter(const ForestStorage& forest_storage, std::ostream& full, std::ostream& included);

} // End namespace ig_simulator
