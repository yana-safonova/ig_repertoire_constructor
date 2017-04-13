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
    TreeCreator tree_creator;

public:
    ForestCreator(AbstractShmCreatorCPtr&& shm_creator,
                  AbstractTreeSizeGeneratorCPtr&& tree_size_generator,
                  double ret_prob):
        tree_creator(std::move(shm_creator), std::move(tree_size_generator), ret_prob)
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
