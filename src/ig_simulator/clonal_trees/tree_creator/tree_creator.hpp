//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "clonal_trees/tree/tree.hpp"
#include "pool_manager.hpp"
#include "base_repertoire/metaroot/metaroot.hpp"
#include "shm_creator.hpp"
#include "tree_size_generator.hpp"
#include "clonal_trees/fast_stop_codon_checker/fast_stop_codon_checker.hpp"

namespace ig_simulator {

class TreeCreator {
protected:
    const AbstractShmCreatorCPtr shm_creator;
    const AbstractTreeSizeGeneratorCPtr tree_size_generator;
    const double ret_prob;

    mutable std::geometric_distribution<size_t> distr_n_children;

private:
    std::string CreateSequence(const std::string& base_seq, const Node::SHM_Vector& shms) const {
        std::string seq = base_seq;

        for(const auto& shm : shms) {
            VERIFY_MSG(seq[std::get<0>(shm)] == std::get<1>(shm),
                       std::string("real seq: ") << seq <<
                       ", position: " << std::get<0>(shm) <<
                       ", expected: " << std::get<1>(shm));
            seq[std::get<0>(shm)] = std::get<2>(shm);
        }
        return seq;
    }

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

    TreeCreator(const vj_finder::VJFinderConfig& vjf_config,
                const ClonalTreeSimulatorParams& config):
        TreeCreator(get_shm_creator(vjf_config, config.shm_creator_params),
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
        nodes.emplace_back();

        std::vector<std::string> sequences;
        sequences.reserve(tree_size);
        sequences.emplace_back(root->Sequence());

        if (not root->IsProductive()) {
            nodes.back().MakeNonProductive();
            return Tree(root, std::move(nodes), std::move(sequences));
        }

        PoolManager pool_manager(ret_prob);

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
                const std::string& base_sequence = sequences[parent_ind];
                Node::SHM_Vector shm_vector { shm_creator->GenerateSHM_Vector(base_sequence)};
                std::string sequence = CreateSequence(base_sequence, shm_vector);

                nodes.emplace_back(parent_ind, std::move(shm_vector));
                sequences.emplace_back(std::move(sequence));

                if (FastStopCodonChecker::HasStopCodon(sequences.back(), root->CDRLabeling())) {
                    nodes.back().MakeNonProductive();
                    pool_manager.Erase(pool_manager.MaxIndex() - n_children + i);
                }
            }
            if (pool_manager.Size() == 0) { break; } // All leafs are non-productive
        }
        if (pool_manager.Size() != 0) { // Only when leafs all non-productive VERIFY should not be checked.
            VERIFY(nodes.size() == tree_size);
        }
        return Tree(root, std::move(nodes), std::move(sequences));
    }
};

} // End namespace ig_simulator
