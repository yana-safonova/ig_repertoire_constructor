//
// Created by Andrew Bzikadze on 4/14/17.
//

#include "exporters.hpp"

namespace ig_simulator {

void TreeExporter(const Tree& tree, size_t forest_ind, size_t tree_ind,
                         std::ostream& full, std::ostream& included)
{
    const auto sequences = tree.Sequences();
    for (size_t i = 0; i < sequences.size(); ++i) {
        std::stringstream id_ss;
        id_ss << ">forest_" << forest_ind << "_tree_" << tree_ind << "_antibody_" << i;
        std::string id { id_ss.str() };
        full << id << '\n' << sequences[i] << '\n';
        if (tree.IsNodeIncluded(i)) {
            included << id << '\n' << sequences[i] << '\n';
        }
    }
}

void ForestExporter(const Forest& forest, size_t forest_ind, std::ostream& full, std::ostream& included) {
    for (size_t i = 0; i < forest.Trees().size(); ++i) {
        TreeExporter(forest.Trees()[i], forest_ind, i, full, included);
    }
}

void ForestStorageExporter(const ForestStorage& forest_storage, std::ostream& full, std::ostream& included) {
    for (size_t i = 0; i < forest_storage.size(); ++i) {
        ForestExporter(forest_storage[i], i, full, included);
    }
}

void EdgeListsExporters(const ForestStorage& forest_storage, const IgSimulatorConfig::IOParams::OutputParams& config) {
    std::string path = path::append_path(config.output_dir, config.trees_dir);
    path::make_dir(path);
    for (size_t i = 0; i < forest_storage.size(); ++i) {
        for (size_t j = 0; j < forest_storage[i].Size(); ++j) {
            std::stringstream filename;
            filename << "forest_" << i << "_tree_" << j << ".dot";
            std::string full_filename = path::append_path(path, filename.str());
            std::ofstream out;
            out.open(full_filename);
            out << forest_storage[i].Trees()[j];
            out.close();
        }
    }
}

} // End namespace ig_simulator
