//
// Created by Andrew Bzikadze on 4/14/17.
//

#pragma once

#include "forest_creator.hpp"
#include "base_repertoire/base_repertoire.hpp"

namespace ig_simulator {

class ForestStorageCreator {
private:
    ForestCreator forest_creator;

public:
    ForestStorageCreator(AbstractShmCreatorCPtr&& shm_creator,
                         AbstractTreeSizeGeneratorCPtr&& tree_size_generator,
                         double ret_prob):
        forest_creator(std::move(shm_creator), std::move(tree_size_generator), ret_prob)
    { }

    ForestStorageCreator(const ForestStorageCreator&) = delete;
    ForestStorageCreator(ForestStorageCreator&&) = delete;
    ForestStorageCreator& operator=(const ForestStorageCreator&) = delete;
    ForestStorageCreator& operator=(ForestStorageCreator&&) = delete;

    template<class PoolManager>
    ForestStorage GenerateForest(const BaseRepertoire& repertoire) const {
        ForestStorage storage;
        for(const auto& cluster : repertoire) {
            storage.emplace_back(forest_creator.GenerateForest<PoolManager>(cluster));
        }
        return storage;
    }
};

} // End namespace ig_simulator
