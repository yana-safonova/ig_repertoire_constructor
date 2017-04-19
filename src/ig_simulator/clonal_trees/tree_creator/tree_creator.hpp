//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "clonal_trees/tree/tree.hpp"
#include "pool_manager.hpp"
#include "base_repertoire/metaroot/metaroot.hpp"
#include "shm_creator.hpp"
#include "tree_size_generator.hpp"

namespace ig_simulator {

class TreeCreator {
protected:
    AbstractShmCreatorCPtr shm_creator;
    AbstractTreeSizeGeneratorCPtr tree_size_generator;
    double ret_prob;
    mutable std::geometric_distribution<size_t> distr_n_children;

public:
    TreeCreator(AbstractShmCreatorCPtr&& shm_creator,
                AbstractTreeSizeGeneratorCPtr&& tree_size_generator,
                double ret_prob,
                double lambda_distr_n_children):
        shm_creator(std::move(shm_creator)),
        tree_size_generator(std::move(tree_size_generator)),
        ret_prob(check_numeric_positive(ret_prob)),
        distr_n_children(check_numeric_positive(lambda_distr_n_children))
    { }

    TreeCreator(const ClonalTreeSimulatorParams& config):
        TreeCreator(get_shm_creator(config.shm_creator_params),
                    get_tree_size_generator(config.tree_size_generator_params),
                    config.prob_ret_to_pool,
                    config.lambda_distr_n_children)
    { }

    TreeCreator(const TreeCreator&) = delete;
    TreeCreator(TreeCreator&&) = delete;
    TreeCreator& operator=(const TreeCreator&) = delete;
    TreeCreator& operator=(TreeCreator&&) = delete;

    template<class PoolManager>
    Tree GenerateTree(const AbstractMetaroot* const root) const {
        static_assert(std::is_base_of<AbstractPoolManager, PoolManager>::value,
                      "Pool Manager should be derived from @class AbstractPoolManager");

        size_t tree_size = tree_size_generator->Generate();
        std::vector<Node> nodes;
        nodes.reserve(tree_size);

        PoolManager pool_manager(ret_prob);
        nodes.emplace_back();

        while(nodes.size() < tree_size) {
            size_t n_children = distr_n_children(MTSingleton::GetInstance()) + 1;
            n_children = std::min(n_children, tree_size - nodes.size());
            size_t parent_ind;
            bool stay;
            std::tie(parent_ind, stay) = pool_manager.GetIndex(n_children);

            if (not stay) {
                nodes[parent_ind].Exclude();
            }

            for (size_t i = 0; i < n_children; ++i) {
                Node::SHM_Vector shm_vector { shm_creator->GenerateSHM_Vector(root->Length())};
                nodes.emplace_back(parent_ind, std::move(shm_vector));
            }
        }
        VERIFY(nodes.size() == tree_size);
        return Tree(root, std::move(nodes));
    }
};

} // End namespace ig_simulator
