//
// Created by Andrew Bzikadze on 4/14/17.
//

#pragma once

#include "tree_creator.hpp"
#include "clonal_trees/forest/forest.hpp"
#include "base_repertoire/metaroot_cluster/metaroot_cluster.hpp"

namespace ig_simulator {

class ForestCreator {
private:
    const TreeCreator tree_creator;

public:
    ForestCreator(const vj_finder::VJFinderConfig& vjf_config,
                  const ClonalTreeSimulatorParams& config):
        tree_creator(vjf_config, config)
    { }

    ForestCreator(const ForestCreator&) = delete;
    ForestCreator(ForestCreator&&) = delete;
    ForestCreator& operator=(const ForestCreator&) = delete;
    ForestCreator& operator=(ForestCreator&&) = delete;

    template<class PoolManager>
    Forest GenerateForest(const MetarootCluster& root) const {
        std::vector<Tree> trees;
        trees.reserve(root.Multiplicity());
        for(size_t i = 0; i < root.Multiplicity(); ++i) {
            Tree tree { tree_creator.GenerateTree<PoolManager>(root.MetarootPtr().get()) };
            trees.emplace_back(std::move(tree));
        }
        return Forest(&root, std::move(trees));
    }
};

} // End namespace ig_simulator
