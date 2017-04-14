//
// Created by Andrew Bzikadze on 4/14/17.
//

#include "exporters.hpp"

namespace ig_simulator {

void TreeExport(const Tree& tree, size_t forest_ind, size_t tree_ind, std::ostream& full, std::ostream& included) {
    const auto sequences = tree.Sequences();
    for (size_t i = 0; i < sequences.size(); ++i) {
        std::stringstream id_ss;
        id_ss << "forest_" << forest_ind << "_tree_" << tree_ind << "_antibody_" << i;
        std::string id { id_ss.str() };
        seqan::writeRecord(full, id, sequences[i], seqan::Fasta());
        if (tree.IsNodeIncluded(i)) {
            seqan::writeRecord(included, id, sequences[i], seqan::Fasta());
        }
    }
}

void ForestExporter(const Forest& forest, size_t forest_ind, std::ostream& full, std::ostream& included) {
    for (size_t i = 0; i < forest.Trees().size(); ++i) {
        TreeExport(forest.Trees()[i], forest_ind, i, full, included);
    }
}

void ForestStorageExporter(const ForestStorage& forest_storage, std::ostream& full, std::ostream& included) {
    for (size_t i = 0; i < forest_storage.size(); ++i) {
        ForestExporter(forest_storage[i], i, full, included);
    }
}

} // End namespace ig_simulator
